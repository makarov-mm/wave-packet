// WavePacket — 2D time-dependent Schrödinger equation
// Split-step Fourier method:  psi -> e^{-iV dt/2}  F^-1 e^{-i k^2 dt/2} F  e^{-iV dt/2} psi
// Own iterative radix-2 FFT, multithreaded rows + transpose. hbar = m = 1, dx = 1.
// Phase -> hue, amplitude -> brightness. Draw potential walls with the mouse.
// C++17 / Qt6 Widgets, no other dependencies.

#include <QApplication>
#include <QWidget>
#include <QPainter>
#include <QTimer>
#include <QElapsedTimer>
#include <QMouseEvent>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QLabel>
#include <QGroupBox>
#include <QStyleFactory>

#include <cstdlib>
#include <cmath>
#include <complex>
#include <cstdint>
#include <thread>
#include <vector>

using cfloat = std::complex<float>;
static constexpr float PI = 3.14159265358979f;

static inline float clampf(float v, float lo, float hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

// ----------------------------------------------------------------------- FFT

// In-place iterative radix-2 FFT over contiguous rows.
class Fft
{
public:
    void prepare(int n)
    {
        if (n == m_n) return;
        m_n = n;
        m_rev.resize(n);
        int lg = 0;
        while ((1 << lg) < n) ++lg;
        for (int i = 0; i < n; ++i) {
            int r = 0;
            for (int b = 0; b < lg; ++b)
                if (i & (1 << b)) r |= 1 << (lg - 1 - b);
            m_rev[i] = r;
        }
        m_twF.resize(n / 2);
        m_twI.resize(n / 2);
        for (int i = 0; i < n / 2; ++i) {
            float a = -2.f * PI * i / n;
            m_twF[i] = { std::cos(a),  std::sin(a) };
            m_twI[i] = { std::cos(a), -std::sin(a) };
        }
    }

    // sign: -1 forward, +1 inverse (inverse also scales by 1/n)
    void run(cfloat *a, int sign) const
    {
        const int n = m_n;
        for (int i = 0; i < n; ++i) {
            int j = m_rev[i];
            if (i < j) std::swap(a[i], a[j]);
        }
        const cfloat *tw = sign < 0 ? m_twF.data() : m_twI.data();
        for (int len = 2; len <= n; len <<= 1) {
            const int half = len >> 1;
            const int stride = n / len;
            for (int i = 0; i < n; i += len) {
                for (int k = 0; k < half; ++k) {
                    cfloat w = tw[k * stride];
                    cfloat u = a[i + k];
                    cfloat v = a[i + k + half] * w;
                    a[i + k]        = u + v;
                    a[i + k + half] = u - v;
                }
            }
        }
        if (sign > 0) {
            const float inv = 1.f / n;
            for (int i = 0; i < n; ++i) a[i] *= inv;
        }
    }

private:
    int m_n = 0;
    std::vector<int> m_rev;
    std::vector<cfloat> m_twF, m_twI;
};

// ----------------------------------------------------------------- simulation

enum class Scenario { FreePacket, DoubleSlit, TunnelBarrier, HarmonicTrap, DrawnOnly };

struct SimParams {
    Scenario scenario = Scenario::DoubleSlit;
    int   gridExp     = 9;      // 2^9 = 512
    float dt          = 0.12f;
    int   stepsPerFrame = 3;
    float k0          = 1.1f;   // initial momentum, rad/sample
    float sigma       = 12.f;   // packet width, samples
    float barrier     = 0.65f;  // barrier height (Free-packet kinetic scale: k0^2/2)
    float wallWidth   = 4.f;    // px
    float slitWidth   = 8.f;    // px
    float slitSep     = 26.f;   // px, center-to-center
    float brightness  = 1.6f;
    float gamma       = 0.6f;
    bool  probabilityOnly = false;
    bool  showPotential   = true;
    bool  paused          = false;
};

class QuantumCanvas : public QWidget
{
    Q_OBJECT
public:
    explicit QuantumCanvas(QWidget *parent = nullptr) : QWidget(parent)
    {
        setMinimumSize(560, 560);
        setMouseTracking(false);
        reset();
        connect(&m_timer, &QTimer::timeout, this, &QuantumCanvas::step);
        m_timer.start(16);
    }

    SimParams params;

    void reset()
    {
        m_N = 1 << params.gridExp;
        const size_t total = size_t(m_N) * m_N;
        m_psi.assign(total, cfloat(0.f, 0.f));
        m_scratch.assign(total, cfloat(0.f, 0.f));
        m_userV.assign(total, 0.f);
        m_image = QImage(m_N, m_N, QImage::Format_RGB32);
        m_fft.prepare(m_N);
        rebuildPotential();
        rebuildKineticPhase();
        rebuildAbsorber();
        launchPacket();
    }

    void relaunch() { launchPacket(); }

    void rebuildPotential()
    {
        const int N = m_N;
        m_V.assign(size_t(N) * N, 0.f);
        const float wallV = 4.0f; // tall enough to reflect at these energies
        auto V = [&](int x, int y) -> float & { return m_V[size_t(y) * N + x]; };

        switch (params.scenario) {
        case Scenario::DoubleSlit: {
            const int wx0 = int(N * 0.52f - params.wallWidth * 0.5f);
            const int wx1 = int(N * 0.52f + params.wallWidth * 0.5f);
            const float cy = N * 0.5f;
            const float g0lo = cy - params.slitSep * 0.5f - params.slitWidth * 0.5f;
            const float g0hi = cy - params.slitSep * 0.5f + params.slitWidth * 0.5f;
            const float g1lo = cy + params.slitSep * 0.5f - params.slitWidth * 0.5f;
            const float g1hi = cy + params.slitSep * 0.5f + params.slitWidth * 0.5f;
            for (int y = 0; y < N; ++y) {
                bool gap = (y >= g0lo && y <= g0hi) || (y >= g1lo && y <= g1hi);
                if (gap) continue;
                for (int x = wx0; x <= wx1 && x < N; ++x)
                    if (x >= 0) V(x, y) = wallV;
            }
            break;
        }
        case Scenario::TunnelBarrier: {
            const float E = 0.5f * params.k0 * params.k0;
            const float Vb = params.barrier * (E > 1e-4f ? E / 0.5f : 1.f);
            const int wx0 = int(N * 0.52f - params.wallWidth * 0.5f);
            const int wx1 = int(N * 0.52f + params.wallWidth * 0.5f);
            for (int y = 0; y < N; ++y)
                for (int x = std::max(0, wx0); x <= wx1 && x < N; ++x)
                    V(x, y) = Vb;
            break;
        }
        case Scenario::HarmonicTrap: {
            const float w2 = 8.f / (N * N); // gentle trap
            const float c = N * 0.5f;
            for (int y = 0; y < N; ++y)
                for (int x = 0; x < N; ++x) {
                    float dx = x - c, dy = y - c;
                    V(x, y) = 0.5f * w2 * (dx * dx + dy * dy);
                }
            break;
        }
        case Scenario::FreePacket:
        case Scenario::DrawnOnly:
            break;
        }

        // add user-drawn walls
        const size_t total = size_t(N) * N;
        for (size_t i = 0; i < total; ++i)
            m_V[i] += m_userV[i];
    }

    void rebuildKineticPhase()
    {
        const int N = m_N;
        m_kin.resize(size_t(N) * N);
        std::vector<float> k2(N);
        for (int i = 0; i < N; ++i) {
            int f = i <= N / 2 ? i : i - N;
            float k = 2.f * PI * f / N;
            k2[i] = k * k;
        }
        const float dt = params.dt;
        for (int y = 0; y < N; ++y)
            for (int x = 0; x < N; ++x) {
                float ph = -0.5f * (k2[x] + k2[y]) * dt;
                m_kin[size_t(y) * N + x] = { std::cos(ph), std::sin(ph) };
            }
    }

    void clearDrawn()
    {
        std::fill(m_userV.begin(), m_userV.end(), 0.f);
        rebuildPotential();
    }

    float msPerStep() const { return m_msPerStep; }
    float norm() const { return m_norm; }

signals:
    void statsChanged();

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter p(this);
        p.fillRect(rect(), QColor(6, 6, 10));
        const int side = std::min(width(), height());
        m_view = QRect((width() - side) / 2, (height() - side) / 2, side, side);
        p.setRenderHint(QPainter::SmoothPixmapTransform);
        p.drawImage(m_view, m_image);
    }

    void mousePressEvent(QMouseEvent *e) override { drawAt(e->pos(), e->buttons()); }
    void mouseMoveEvent(QMouseEvent *e) override  { drawAt(e->pos(), e->buttons()); }

private slots:
    void step()
    {
        if (!params.paused) {
            QElapsedTimer t; t.start();
            for (int s = 0; s < params.stepsPerFrame; ++s)
                splitStep();
            m_msPerStep = float(t.nsecsElapsed()) * 1e-6f / params.stepsPerFrame;
        }
        composeImage();
        update();
        if (const char *df = getenv("DUMP_FRAMES")) {
            static int left = atoi(df);
            if (--left <= 0) { m_image.save("dump.png"); QApplication::quit(); }
        }
        if (!m_statClock.isValid() || m_statClock.hasExpired(400)) {
            computeNorm();
            emit statsChanged();
            m_statClock.restart();
        }
    }

private:
    static int threadCount()
    {
        return std::max(2u, std::thread::hardware_concurrency());
    }

    template <class F> void parallelRows(F &&fn) const
    {
        const int T = threadCount();
        std::vector<std::thread> th;
        for (int t = 0; t < T; ++t)
            th.emplace_back([&, t] { for (int y = t; y < m_N; y += T) fn(y); });
        for (auto &x : th) x.join();
    }

    void drawAt(QPoint pos, Qt::MouseButtons b)
    {
        if (!m_view.contains(pos) || !(b & (Qt::LeftButton | Qt::RightButton)))
            return;
        const int N = m_N;
        const float fx = (pos.x() - m_view.x()) / float(m_view.width());
        const float fy = (pos.y() - m_view.y()) / float(m_view.height());
        const int cx = int(fx * N), cy = int(fy * N);
        const int r = std::max(2, N / 64);
        const float v = (b & Qt::LeftButton) ? 4.0f : 0.f;
        for (int y = std::max(0, cy - r); y <= std::min(N - 1, cy + r); ++y)
            for (int x = std::max(0, cx - r); x <= std::min(N - 1, cx + r); ++x) {
                int dx = x - cx, dy = y - cy;
                if (dx * dx + dy * dy <= r * r)
                    m_userV[size_t(y) * N + x] = v;
            }
        rebuildPotential();
    }

    void launchPacket()
    {
        const int N = m_N;
        float x0 = N * 0.28f, y0 = N * 0.5f;
        float kx = params.k0, ky = 0.f;
        if (params.scenario == Scenario::HarmonicTrap) {
            x0 = N * 0.36f; // displaced -> classical oscillation
            kx = 0.f; ky = params.k0 * 0.4f;
        }
        const float s2 = 2.f * params.sigma * params.sigma;
        parallelRows([&](int y) {
            cfloat *row = &m_psi[size_t(y) * N];
            for (int x = 0; x < N; ++x) {
                float dx = x - x0, dy = y - y0;
                float amp = std::exp(-(dx * dx + dy * dy) / s2);
                float ph = kx * x + ky * y;
                row[x] = cfloat(amp * std::cos(ph), amp * std::sin(ph));
            }
        });
        normalize();
    }

    void normalize()
    {
        double sum = 0.0;
        for (const auto &c : m_psi) sum += double(std::norm(c));
        if (sum < 1e-20) return;
        const float inv = float(1.0 / std::sqrt(sum));
        for (auto &c : m_psi) c *= inv;
    }

    void computeNorm()
    {
        double sum = 0.0;
        for (const auto &c : m_psi) sum += double(std::norm(c));
        m_norm = float(sum);
    }

    void rebuildAbsorber()
    {
        const int N = m_N;
        const int W = std::max(12, N / 20);
        m_absorb.resize(N);
        for (int i = 0; i < N; ++i) {
            int d = std::min(i, N - 1 - i);
            if (d >= W) { m_absorb[i] = 1.f; continue; }
            float t = float(d) / W;
            float ramp = 0.5f - 0.5f * std::cos(PI * t);
            m_absorb[i] = 0.965f + 0.035f * ramp; // gentle sponge
        }
    }

    void applyPotentialHalfStep()
    {
        const float h = -0.5f * params.dt;
        parallelRows([&](int y) {
            cfloat *row = &m_psi[size_t(y) * m_N];
            const float *v = &m_V[size_t(y) * m_N];
            const float ay = m_absorb[y];
            for (int x = 0; x < m_N; ++x) {
                float ph = h * v[x];
                row[x] *= cfloat(std::cos(ph), std::sin(ph)) * (ay * m_absorb[x]);
            }
        });
    }

    void transpose(std::vector<cfloat> &src, std::vector<cfloat> &dst) const
    {
        const int N = m_N, B = 32;
        const int T = threadCount();
        std::vector<std::thread> th;
        for (int t = 0; t < T; ++t) {
            th.emplace_back([&, t] {
                for (int by = t * B; by < N; by += T * B)
                    for (int bx = 0; bx < N; bx += B)
                        for (int y = by; y < std::min(N, by + B); ++y)
                            for (int x = bx; x < std::min(N, bx + B); ++x)
                                dst[size_t(x) * N + y] = src[size_t(y) * N + x];
            });
        }
        for (auto &x : th) x.join();
    }

    void fftRows(std::vector<cfloat> &data, int sign)
    {
        parallelRows([&](int y) { m_fft.run(&data[size_t(y) * m_N], sign); });
    }

    void splitStep()
    {
        applyPotentialHalfStep();

        // 2D FFT: rows, transpose, rows, transpose back
        fftRows(m_psi, -1);
        transpose(m_psi, m_scratch);
        fftRows(m_scratch, -1);
        transpose(m_scratch, m_psi); // psi now in k-space, natural layout

        parallelRows([&](int y) {
            cfloat *row = &m_psi[size_t(y) * m_N];
            const cfloat *k = &m_kin[size_t(y) * m_N];
            for (int x = 0; x < m_N; ++x) row[x] *= k[x];
        });

        fftRows(m_psi, +1);
        transpose(m_psi, m_scratch);
        fftRows(m_scratch, +1);
        transpose(m_scratch, m_psi);

        applyPotentialHalfStep();
    }

    void composeImage()
    {
        const int N = m_N;
        const bool prob = params.probabilityOnly;
        const bool showV = params.showPotential;
        const float bright = params.brightness * N * 0.06f;
        const float gamma = params.gamma;
        parallelRows([&](int y) {
            QRgb *out = reinterpret_cast<QRgb *>(m_image.scanLine(y));
            const cfloat *row = &m_psi[size_t(y) * N];
            const float *v = &m_V[size_t(y) * N];
            for (int x = 0; x < N; ++x) {
                float a = std::abs(row[x]);
                float val = clampf(std::pow(a * bright, gamma), 0.f, 1.f);
                int r, g, b;
                if (prob) {
                    // inferno-ish ramp
                    r = int(clampf(val * 2.4f, 0.f, 1.f) * 255);
                    g = int(clampf(val * 1.5f - 0.25f, 0.f, 1.f) * 255);
                    b = int(clampf(0.4f * val + (val > 0.85f ? (val - 0.85f) * 4.f : 0.f), 0.f, 1.f) * 255);
                } else {
                    float hue = (std::atan2(row[x].imag(), row[x].real()) + PI) / (2.f * PI);
                    float h6 = hue * 6.f;
                    int   hi = int(h6) % 6;
                    float f = h6 - int(h6);
                    float p = 0.f, q = 1.f - f, t = f;
                    float rr, gg, bb;
                    switch (hi) {
                    case 0: rr = 1; gg = t; bb = p; break;
                    case 1: rr = q; gg = 1; bb = p; break;
                    case 2: rr = p; gg = 1; bb = t; break;
                    case 3: rr = p; gg = q; bb = 1; break;
                    case 4: rr = t; gg = p; bb = 1; break;
                    default: rr = 1; gg = p; bb = q; break;
                    }
                    r = int(rr * val * 255);
                    g = int(gg * val * 255);
                    b = int(bb * val * 255);
                }
                if (showV && v[x] > 0.01f) {
                    int w = int(clampf(v[x] * 0.25f, 0.08f, 0.55f) * 255);
                    r = std::min(255, r + w);
                    g = std::min(255, g + w / 3);
                    b = std::min(255, b + w / 4);
                }
                out[x] = qRgb(r, g, b);
            }
        });
    }

    int m_N = 512;
    std::vector<cfloat> m_psi, m_scratch, m_kin;
    std::vector<float> m_V, m_userV, m_absorb;
    Fft m_fft;
    QImage m_image;
    QRect m_view;
    QTimer m_timer;
    QElapsedTimer m_statClock;
    float m_msPerStep = 0.f, m_norm = 1.f;
};

// -------------------------------------------------------------------- window

class MainWindow : public QWidget
{
    Q_OBJECT
public:
    MainWindow()
    {
        setWindowTitle("WavePacket · 2D Schrodinger solver");
        auto *canvas = new QuantumCanvas;

        auto *scenario = new QComboBox;
        scenario->addItems({ "Free packet", "Double slit", "Tunnel barrier",
                             "Harmonic trap", "Empty (draw walls)" });
        scenario->setCurrentIndex(1);

        auto *grid = new QComboBox;
        grid->addItems({ "256 x 256", "512 x 512", "1024 x 1024" });
        grid->setCurrentIndex(1);

        auto *dt = new QDoubleSpinBox; dt->setRange(0.01, 0.5); dt->setDecimals(3);
        dt->setSingleStep(0.01); dt->setValue(canvas->params.dt);
        auto *steps = new QSpinBox; steps->setRange(1, 20); steps->setValue(canvas->params.stepsPerFrame);
        auto *k0 = new QDoubleSpinBox; k0->setRange(0.0, 2.2); k0->setSingleStep(0.05);
        k0->setValue(canvas->params.k0);
        auto *sigma = new QDoubleSpinBox; sigma->setRange(3.0, 60.0); sigma->setSingleStep(1.0);
        sigma->setValue(canvas->params.sigma);
        auto *barrier = new QDoubleSpinBox; barrier->setRange(0.0, 3.0); barrier->setSingleStep(0.05);
        barrier->setValue(canvas->params.barrier);
        auto *slitW = new QDoubleSpinBox; slitW->setRange(2.0, 40.0); slitW->setValue(canvas->params.slitWidth);
        auto *slitS = new QDoubleSpinBox; slitS->setRange(6.0, 120.0); slitS->setValue(canvas->params.slitSep);
        auto *bright = new QDoubleSpinBox; bright->setRange(0.2, 8.0); bright->setSingleStep(0.1);
        bright->setValue(canvas->params.brightness);

        auto *probOnly = new QCheckBox("Probability only");
        auto *showV = new QCheckBox("Show potential");
        showV->setChecked(true);

        auto *relaunch = new QPushButton("Relaunch packet");
        auto *clearWalls = new QPushButton("Clear drawn walls");
        auto *pause = new QPushButton("Pause");
        pause->setCheckable(true);

        auto *hint = new QLabel("LMB: draw wall\nRMB: erase");
        hint->setStyleSheet("color:#889");
        auto *stats = new QLabel;
        stats->setStyleSheet("color:#8fa;font-family:monospace");

        auto *form = new QFormLayout;
        form->addRow("Scenario", scenario);
        form->addRow("Grid", grid);
        form->addRow("dt", dt);
        form->addRow("Steps/frame", steps);
        form->addRow("Momentum k0", k0);
        form->addRow("Packet sigma", sigma);
        form->addRow("Barrier V/E", barrier);
        form->addRow("Slit width", slitW);
        form->addRow("Slit separation", slitS);
        form->addRow("Brightness", bright);
        form->addRow(probOnly);
        form->addRow(showV);
        form->addRow(relaunch);
        form->addRow(clearWalls);
        form->addRow(pause);
        form->addRow(hint);
        form->addRow(stats);

        auto *group = new QGroupBox("Parameters");
        group->setLayout(form);
        group->setFixedWidth(280);

        auto *layout = new QHBoxLayout(this);
        layout->addWidget(group);
        layout->addWidget(canvas, 1);

        connect(scenario, &QComboBox::currentIndexChanged, this, [canvas](int i) {
            canvas->params.scenario = Scenario(i);
            canvas->rebuildPotential();
            canvas->relaunch();
        });
        connect(grid, &QComboBox::currentIndexChanged, this, [canvas](int i) {
            canvas->params.gridExp = 8 + i;
            canvas->reset();
        });
        connect(dt, &QDoubleSpinBox::valueChanged, this, [canvas](double v) {
            canvas->params.dt = float(v);
            canvas->rebuildKineticPhase();
        });
        connect(steps, &QSpinBox::valueChanged, this, [canvas](int v) { canvas->params.stepsPerFrame = v; });
        connect(k0, &QDoubleSpinBox::valueChanged, this, [canvas](double v) {
            canvas->params.k0 = float(v);
            canvas->rebuildPotential(); // barrier height tracks E
        });
        connect(sigma, &QDoubleSpinBox::valueChanged, this, [canvas](double v) { canvas->params.sigma = float(v); });
        connect(barrier, &QDoubleSpinBox::valueChanged, this, [canvas](double v) {
            canvas->params.barrier = float(v);
            canvas->rebuildPotential();
        });
        connect(slitW, &QDoubleSpinBox::valueChanged, this, [canvas](double v) {
            canvas->params.slitWidth = float(v);
            canvas->rebuildPotential();
        });
        connect(slitS, &QDoubleSpinBox::valueChanged, this, [canvas](double v) {
            canvas->params.slitSep = float(v);
            canvas->rebuildPotential();
        });
        connect(bright, &QDoubleSpinBox::valueChanged, this, [canvas](double v) { canvas->params.brightness = float(v); });
        connect(probOnly, &QCheckBox::toggled, this, [canvas](bool b) { canvas->params.probabilityOnly = b; });
        connect(showV, &QCheckBox::toggled, this, [canvas](bool b) { canvas->params.showPotential = b; });
        connect(relaunch, &QPushButton::clicked, this, [canvas] { canvas->relaunch(); });
        connect(clearWalls, &QPushButton::clicked, this, [canvas] { canvas->clearDrawn(); });
        connect(pause, &QPushButton::toggled, this, [canvas, pause](bool b) {
            canvas->params.paused = b;
            pause->setText(b ? "Resume" : "Pause");
        });
        connect(canvas, &QuantumCanvas::statsChanged, this, [canvas, stats] {
            stats->setText(QString("%1 ms/step\nnorm %2")
                           .arg(canvas->msPerStep(), 0, 'f', 2)
                           .arg(canvas->norm(), 0, 'f', 4));
        });

        resize(1180, 900);
    }
};

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    app.setStyle(QStyleFactory::create("Fusion"));
    QPalette pal;
    pal.setColor(QPalette::Window, QColor(37, 37, 42));
    pal.setColor(QPalette::WindowText, QColor(220, 220, 224));
    pal.setColor(QPalette::Base, QColor(28, 28, 32));
    pal.setColor(QPalette::Text, QColor(220, 220, 224));
    pal.setColor(QPalette::Button, QColor(48, 48, 54));
    pal.setColor(QPalette::ButtonText, QColor(220, 220, 224));
    pal.setColor(QPalette::Highlight, QColor(70, 120, 200));
    app.setPalette(pal);

    MainWindow w;
    w.show();
    return app.exec();
}

#include "main.moc"

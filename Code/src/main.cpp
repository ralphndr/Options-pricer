#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <memory>
#include <chrono>
// formatting
#include <iomanip>
// Project headers
#include "MarketData.hpp"
#include "VanillaOption.hpp"
#include "ExoticOption.hpp"
#include "ExplicitSolver.hpp"
#include "ImplicitSolver.hpp"
#include "CrankNicolsonSolver.hpp"
#include "AmericanOption.hpp"
#include "AmericanSolver.hpp"
#include "MonteCarlo.hpp"
 

struct Args {
    std::string optionClass = "european"; // european | american
    std::string optionType  = "call";     // call | put
    std::string solver      = "cn";       // explicit | implicit | cn | american | mc
    int M = 200;
    int N = 2000;
    double Smax = 150.0;
    double spot = 100.0;    // spot for Monte Carlo
    double strike = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;
    bool plot = true;
    bool viewer = false;   // viewer mode
    bool spin = false;     // rotate view
    double clampK = -1.0;  // clamp z/color range to [0:K]
    std::string out = "output_grid.csv";
    std::string earlyBoundaryOut = "";
    bool overlayEarlyBoundary = false;
    std::string barrierType = "";
    double barrierLevel = -1.0;
    double digitalCash = 1.0;
    int mcPaths = 50000;
    int mcSteps = 252;
    unsigned int mcSeed = 0; // 0 => random_device
    bool mcAntithetic = true;
    double mcJumpIntensity = 0.0;
    double mcJumpMean = 0.0;
    double mcJumpVol = 0.0;
    std::string lookbackType = "";
    double chooseTime = -1.0;
    bool printPlotData = false; // print space-delimited plot data to stdout
    int printCount = 1000;      // number of lines to print
    // Plot export options
    std::string plotFormat = ""; // "png" or "pdf"; empty => GUI viewer
    std::string plotOut = "";    // output image path; default based on format
    int plotWidth = 1200;         // pixels for PNG; logical size for viewer
    int plotHeight = 800;         // pixels for PNG; logical size for viewer
    int plotDpi = 150;            // PNG only; ignored for PDF
};

static void print_usage() {
    std::cout << "Usage: option_viz [flags]\n"
              << "  --option european|american\n"
              << "  --type call|put\n"
              << "  --solver explicit|implicit|cn|american|mc\n"
              << "  --M <int> --N <int> --Smax <float>\n"
              << "  --spot <float> (for Monte Carlo)\n"
              << "  --strike <float> --T <float> --r <float> --sigma <float>\n"
              << "  --out <path>\n"
              << "  --plot [on|off]\n"
              << "  --viewer (open plot directly with current args)\n"
              << "  --spin [on|off] (rotate the 3D view)\n"
              << "  --clampK <float> (fix z/color to [0:K])\n"
              << "  --earlyBoundaryOut <path> (export S*(tau), V*(tau) for American)\n"
              << "  --overlayEarlyBoundary [on|off] (overlay S*(tau) curve in viewer)\n"
              << "  --barrierType downout|upout (for barrier options)\n"
              << "  --barrierLevel <float> (barrier level B)\n"
              << "  --digitalCash <float> (cash payout for digital option)\n"
              << "  --mcPaths <int> (Monte Carlo paths)\n"
              << "  --mcSteps <int> (timesteps per path)\n"
              << "  --mcSeed <uint> (0 => random)\n"
              << "  --mcAntithetic [on|off] (antithetic variates)\n"
              << "  --mcJumpIntensity <float> (lambda; 0 disables jumps)\n"
              << "  --mcJumpMean <float> (jump log-mean mu_J)\n"
              << "  --mcJumpVol <float> (jump log-vol sigma_J)\n"
              << "  --lookbackType min|max (for lookback options)\n"
              << "  --chooseTime <float> (choose time as fraction of T, e.g., 0.5)\n"
              << "  --printPlotData [on|off] (print space-delimited t S V lines)\n"
              << "  --printCount <int> (how many lines to print when enabled)\n"
              << "\nPlot export (gnuplot)\n"
              << "  --plotFormat png|pdf (save image instead of opening viewer)\n"
              << "  --plotOut <path> (output image path)\n"
              << "  --plotWidth <int> --plotHeight <int> (image size; PNG pixels, PDF in inches via 96dpi)\n"
              << "  --plotDpi <int> (PNG only; default 150)\n";
}

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto next = [&](double &dest){ if (i+1 < argc) dest = std::stod(argv[++i]); };
        auto nexti = [&](int &dest){ if (i+1 < argc) dest = std::stoi(argv[++i]); };
        auto nexts = [&](std::string &dest){ if (i+1 < argc) dest = argv[++i]; };
        if (k == "--option") nexts(a.optionClass);
        else if (k == "--type") nexts(a.optionType);
        else if (k == "--solver") nexts(a.solver);
        else if (k == "--M") nexti(a.M);
        else if (k == "--N") nexti(a.N);
        else if (k == "--Smax") next(a.Smax);
        else if (k == "--spot") next(a.spot);
        else if (k == "--strike") next(a.strike);
        else if (k == "--T") next(a.T);
        else if (k == "--r") next(a.r);
        else if (k == "--sigma") next(a.sigma);
        else if (k == "--out") nexts(a.out);
        else if (k == "--earlyBoundaryOut") nexts(a.earlyBoundaryOut);
        else if (k == "--overlayEarlyBoundary") {
            std::string v = (i+1 < argc) ? argv[i+1] : "off";
            if (v == "on" || v == "off") { a.overlayEarlyBoundary = (v == "on"); ++i; }
        }
        else if (k == "--barrierType") nexts(a.barrierType);
        else if (k == "--barrierLevel") next(a.barrierLevel);
        else if (k == "--digitalCash") next(a.digitalCash);
        else if (k == "--mcPaths") nexti(a.mcPaths);
        else if (k == "--mcSteps") nexti(a.mcSteps);
        else if (k == "--mcSeed") { double tmp; next(tmp); a.mcSeed = static_cast<unsigned int>(tmp); }
        else if (k == "--mcAntithetic") {
            std::string v = (i+1 < argc) ? argv[i+1] : "on";
            if (v == "on" || v == "off") { a.mcAntithetic = (v == "on"); ++i; }
        }
        else if (k == "--mcJumpIntensity") next(a.mcJumpIntensity);
        else if (k == "--mcJumpMean") next(a.mcJumpMean);
        else if (k == "--mcJumpVol") next(a.mcJumpVol);
        else if (k == "--lookbackType") nexts(a.lookbackType);
        else if (k == "--chooseTime") next(a.chooseTime);
        else if (k == "--printPlotData") {
            std::string v = (i+1 < argc) ? argv[i+1] : "off";
            if (v == "on" || v == "off") { a.printPlotData = (v == "on"); ++i; }
        }
        else if (k == "--printCount") nexti(a.printCount);
        else if (k == "--plot") {
            std::string v = (i+1 < argc) ? argv[i+1] : "on";
            if (v == "on" || v == "off") { a.plot = (v == "on"); ++i; }
            else { a.plot = true; }
        }
        else if (k == "--viewer") {
            a.viewer = true;
        }
        else if (k == "--plotFormat") { nexts(a.plotFormat); }
        else if (k == "--plotOut") { nexts(a.plotOut); }
        else if (k == "--plotWidth") { nexti(a.plotWidth); }
        else if (k == "--plotHeight") { nexti(a.plotHeight); }
        else if (k == "--plotDpi") { nexti(a.plotDpi); }
        else if (k == "--spin") {
            std::string v = (i+1 < argc) ? argv[i+1] : "off";
            if (v == "on" || v == "off") { a.spin = (v == "on"); ++i; }
        }
        else if (k == "--clampK") { next(a.clampK); }
        else if (k == "-h" || k == "--help") { print_usage(); std::exit(0); }
    }
    return a;
}

// Early exercise boundary detection for American Put
static bool write_early_boundary_csv(const Grid& grid, const Option& option, double T, const std::string& path) {
    if (path.empty()) return false;
    const AmericanOption* am = dynamic_cast<const AmericanOption*>(&option);
    if (!am || am->getType() != AmericanOptionType::Put) return false;
    std::ofstream out(path);
    if (!out) return false;
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();
    double dt = grid.getdt();
    out << "tau,Sstar,Vstar\n";
    for (int j = 0; j <= N; ++j) {
        double tau = T - j * dt;
        int ib = -1;
        for (int i = 0; i <= M; ++i) {
            double S = i * dS;
            double diff = grid.get(i, j) - option.payoff(S);
            if (diff > 1e-8) { ib = i; break; }
        }
        if (ib <= 0) {
            out << tau << "," << 0.0 << "," << grid.get(0, j) << "\n";
        } else {
            double S_lo = (ib-1) * dS;
            double S_hi = ib * dS;
            double d_lo = grid.get(ib-1, j) - option.payoff(S_lo);
            double d_hi = grid.get(ib,   j) - option.payoff(S_hi);
            double w = (0.0 - d_lo) / (d_hi - d_lo + 1e-16);
            if (w < 0.0) w = 0.0; if (w > 1.0) w = 1.0;
            double Sstar = S_lo + w * (S_hi - S_lo);
            double Vstar = option.payoff(Sstar);
            out << tau << "," << Sstar << "," << Vstar << "\n";
        }
    }
    return true;
}

static void write_grid_csv(const Grid& grid, double T, const std::string& path) {
    std::ofstream out(path);
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();
    double dt = grid.getdt();
    out << "S,tau,V\n"; // tau = time until maturity
    int stepS = std::max(1, M / 200); // thin for plotting speed
    int stepT = std::max(1, N / 200);
    for (int j = 0; j <= N; j += stepT) {
        double tau = T - j * dt;
        for (int i = 0; i <= M; i += stepS) {
            double S = i * dS;
            out << S << "," << tau << "," << grid.get(i, j) << "\n";
        }
    }
}

// Print space-delimited plot data in the order: t S V (no header)
// t = j*dt (elapsed time), S = i*dS, V = grid(i,j)
static void print_plot_data_stdout(const Grid& grid, int maxLines) {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();
    double dt = grid.getdt();
    int remaining = maxLines;
    for (int j = 0; j <= N && remaining > 0; ++j) {
        double t = j * dt;
        for (int i = 0; i <= M && remaining > 0; ++i) {
            double S = i * dS;
            double V = grid.get(i, j);
            // t and S fixed with 4 decimals; V scientific with 6 decimals
            std::cout << std::fixed << std::setprecision(4) << t << " "
                      << std::fixed << std::setprecision(4) << S << " "
                      << std::scientific << std::setprecision(6) << V << "\n";
            --remaining;
        }
    }
}

// 3D surface using gnuplot
static void plot_surface_gnuplot(const Grid& grid, double T, const std::string& title,
                                 bool spin, double clampK = -1.0,
                                 const std::string& overlayPath = "",
                                 int downS = 200, int downT = 200,
                                 const std::string& termFormat = "",
                                 const std::string& outPath = "",
                                 int width = 1200, int height = 800,
                                 int dpi = 150) {
    (void)overlayPath; // overlay disabled to avoid extra lines on the surface
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();
    double dt = grid.getdt();

    int stepS = std::max(1, M / std::max(1, downS));
    int stepT = std::max(1, N / std::max(1, downT));

    std::string dataPath = "/tmp/option_surface.dat";
    std::ofstream ofs(dataPath);
    if (!ofs) {
        throw std::runtime_error("Failed to open temp data file for gnuplot");
    }
    // Write scattered points: S tau V (matching Python meshgrid order)
    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    for (int j = 0; j <= N; j += stepT) {
        double tau = T - j * dt;
        for (int i = 0; i <= M; i += stepS) {
            double S = i * dS;
            double v = grid.get(i, j);
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
            ofs << S << " " << tau << " " << v << "\n";
        }
        ofs << "\n"; // blank line to separate rows
    }
    ofs.close();

    // Launch gnuplot and send commands
    FILE* gp = popen("gnuplot -persist", "w");
    if (!gp) {
        throw std::runtime_error("Failed to start gnuplot.");
    }
    // Terminal and canvas
    if (termFormat == "png") {
        // pngcairo uses pixel size; DPI implicit in pixel count
        fprintf(gp, "set term pngcairo size %d,%d enhanced\n", width, height);
        std::string path = outPath.empty() ? std::string("option_surface.png") : outPath;
        fprintf(gp, "set output '%s'\n", path.c_str());
    } else if (termFormat == "pdf") {
        // pdfcairo uses physical size; convert pixels to inches via 96 dpi baseline
        double win = std::max(1, width) / 96.0;
        double hin = std::max(1, height) / 96.0;
        fprintf(gp, "set term pdfcairo size %0.3fin,%0.3fin enhanced fontscale 0.9\n", win, hin);
        std::string path = outPath.empty() ? std::string("option_surface.pdf") : outPath;
        fprintf(gp, "set output '%s'\n", path.c_str());
    } else {
        fprintf(gp, "set term qt size %d,%d\n", width, height);
    }
    fprintf(gp, "set view 60,40\n");
    fprintf(gp, "set ticslevel 0\n");
    // Axes ranges and labels
    fprintf(gp, "set xrange [0:200]\n");
    fprintf(gp, "set yrange [0:%g]\n", T + 1e-12);
    fprintf(gp, "set xlabel 'Spot price' offset 0,-1,0\n");
    fprintf(gp, "set ylabel 'Time (in years)' offset 2,0,0\n");
    fprintf(gp, "set zlabel 'Option value' offset 1,0,0\n");
    if (clampK > 0) {
        fprintf(gp, "unset autoscale z\n");
        fprintf(gp, "set zrange [0:%g]\n", clampK + 1e-12);
        fprintf(gp, "set ztics 0,%g,%g\n", clampK/5.0, clampK);
    }
    fprintf(gp, "unset key\n");
    fprintf(gp, "set border lc rgb '#666666' lw 1\n");
    fprintf(gp, "set grid back lc rgb '#bbbbbb' lw 0.5\n");
    fprintf(gp, "set format cb '%%.2f'\n");
    // Color palette (viridis-like)
    fprintf(gp, "set palette defined (\\\n");
    fprintf(gp, " 0 0.267 0.004 0.329,\\\n");
    fprintf(gp, " 0.13 0.283 0.141 0.458,\\\n");
    fprintf(gp, " 0.25 0.254 0.265 0.530,\\\n");
    fprintf(gp, " 0.38 0.206 0.372 0.553,\\\n");
    fprintf(gp, " 0.50 0.163 0.471 0.558,\\\n");
    fprintf(gp, " 0.63 0.127 0.566 0.551,\\\n");
    fprintf(gp, " 0.75 0.135 0.659 0.518,\\\n");
    fprintf(gp, " 0.88 0.267 0.749 0.441,\\\n");
    fprintf(gp, " 1.0 0.993 0.906 0.144)\n");
    fprintf(gp, "set colorbox vertical user origin 0.90,0.20 size 0.03,0.60\n");
    fprintf(gp, "unset autoscale cb\n");
    if (std::isfinite(vmin) && std::isfinite(vmax) && vmax > vmin) {
        fprintf(gp, "set cbrange [%g:%g]\n", vmin, vmax);
        
        // Generate explicit tick labels
        int nticks = 6; // 0, 20, 40, 60, 80, 100 style
        fprintf(gp, "set cbtics (");
        for(int i = 0; i < nticks; i++) {
            double tick_val = vmin + i * (vmax - vmin) / (nticks - 1);
            fprintf(gp, "'%.1f' %g", tick_val, tick_val);
            if(i < nticks - 1) fprintf(gp, ", ");
        }
        fprintf(gp, ")\n");
    }
    fprintf(gp, "set cblabel 'Option value'\n");
    // Surface rendering
    int gx = std::max(40, std::min(160, (M/stepS)+1));
    int gy = std::max(40, std::min(160, (N/stepT)+1));
    fprintf(gp, "set hidden3d\n");
    fprintf(gp, "set pm3d depthorder corners2color mean at s\n");
    // dgrid3d caused errors on some gnuplot builds; omit for broader compatibility
    // Title
    std::string safeTitle = title;
    for (char& c : safeTitle) { if (c=='\'') c=' '; }
    fprintf(gp, "set title '%s' font ',12'\n", safeTitle.c_str());
    // Always render only the surface; omit any overlay line to keep the view clean
    fprintf(gp, "splot '%s' using 1:2:3 with pm3d\n", dataPath.c_str());
    if (spin && termFormat.empty()) {
        // Only spin in interactive viewer
        fprintf(gp, "do for [a=0:360:3] { set view 60,a; replot; pause 0.02 }\n");
    }
    fflush(gp);
    if (!termFormat.empty()) {
        // Ensure file output is finalized
        fprintf(gp, "unset output\n");
        fflush(gp);
    }
    // Do not close gp immediately to keep window open in viewer mode
}


int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);

    // Build market and option
    MarketData market(args.r, args.sigma);

    // Create appropriate option and solver
    std::unique_ptr<Option> opt;
    std::unique_ptr<PDESolver> solver;

    auto to_lower = [](std::string s){ for (auto &c: s) c = std::tolower(c); return s; };
    std::string optClass = to_lower(args.optionClass);
    std::string optType  = to_lower(args.optionType);
    std::string solverNm = to_lower(args.solver);

    // Monte Carlo branch
    if (solverNm == "mc") {
        if (optClass != "european" && optClass != "asian" && optClass != "barrier" && 
            optClass != "lookback" && optClass != "chooser") {
            std::cerr << "Monte Carlo supports option=european|asian|barrier|lookback|chooser\n";
            return 1;
        }
        if (optClass == "barrier" && (args.barrierType.empty() || args.barrierLevel <= 0)) {
            std::cerr << "Barrier MC requires --barrierType and --barrierLevel\n";
            return 1;
        }
        if (optClass == "lookback" && args.lookbackType.empty()) {
            std::cerr << "Lookback MC requires --lookbackType min|max\n";
            return 1;
        }
        OptionType t = (optType == "put") ? OptionType::Put : OptionType::Call;
        MCParams mp;
        mp.paths = args.mcPaths;
        mp.steps = args.mcSteps;
        mp.seed = args.mcSeed;
        mp.antithetic = args.mcAntithetic;
        mp.jumpIntensity = args.mcJumpIntensity;
        mp.jumpMean = args.mcJumpMean;
        mp.jumpVol = args.mcJumpVol;
        mp.lookbackType = args.lookbackType;
        mp.chooseTime = args.chooseTime;

        auto t_start = std::chrono::high_resolution_clock::now();
        MCResult res = mc_price(optClass, t, args.spot, args.strike, args.T, market, mp, args.barrierType, args.barrierLevel);
        auto t_end = std::chrono::high_resolution_clock::now();
        auto elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        std::cout << "MC Price: " << res.price << " (stderr: " << res.stderr << ")\n";
        std::cout << "Time: " << std::fixed << std::setprecision(2) << elapsed_ms << " ms\n";
        return 0;
    }

    if (optClass == "european") {
        OptionType t = (optType == "put") ? OptionType::Put : OptionType::Call;
        opt = std::make_unique<EuropeanOption>(args.strike, args.T, t);
        if (solverNm == "explicit") solver = std::make_unique<ExplicitSolver>(*opt, market, args.M, args.N, args.Smax);
        else if (solverNm == "implicit") solver = std::make_unique<ImplicitSolver>(*opt, market, args.M, args.N, args.Smax);
        else solver = std::make_unique<CrankNicolsonSolver>(*opt, market, args.M, args.N, args.Smax);
    } else if (optClass == "american") {
        AmericanOptionType t = (optType == "put") ? AmericanOptionType::Put : AmericanOptionType::Call;
        auto am = std::make_unique<AmericanOption>(args.strike, args.T, t);
        opt = std::move(am);
        // Force American solver for American option
        solver = std::make_unique<AmericanSolver>(static_cast<const AmericanOption&>(*opt), market, args.M, args.N, args.Smax);
    } else if (optClass == "barrier") {
        OptionType t = (optType == "put") ? OptionType::Put : OptionType::Call;
        if (args.barrierType.empty() || args.barrierLevel <= 0) {
            std::cerr << "Barrier options require --barrierType and --barrierLevel\n";
            return 1;
        }
        // Knock-in via parity: V_in = V_euro - V_out
        bool isIn = (args.barrierType == "downin" || args.barrierType == "upin");
        if (!isIn) {
            opt = std::unique_ptr<Option>(new BarrierOption(args.strike, args.T, t, args.barrierType, args.barrierLevel));
            if (solverNm == "explicit") solver = std::make_unique<ExplicitSolver>(*opt, market, args.M, args.N, args.Smax);
            else if (solverNm == "implicit") solver = std::make_unique<ImplicitSolver>(*opt, market, args.M, args.N, args.Smax);
            else solver = std::make_unique<CrankNicolsonSolver>(*opt, market, args.M, args.N, args.Smax);
        } else {
            // Build out option
            std::string outType = (args.barrierType == "downin") ? "downout" : "upout";
            std::unique_ptr<Option> outOpt(new BarrierOption(args.strike, args.T, t, outType, args.barrierLevel));
            std::unique_ptr<PDESolver> outSol;
            if (solverNm == "explicit") outSol = std::make_unique<ExplicitSolver>(*outOpt, market, args.M, args.N, args.Smax);
            else if (solverNm == "implicit") outSol = std::make_unique<ImplicitSolver>(*outOpt, market, args.M, args.N, args.Smax);
            else outSol = std::make_unique<CrankNicolsonSolver>(*outOpt, market, args.M, args.N, args.Smax);

            // Build vanilla European
            std::unique_ptr<Option> euro(new EuropeanOption(args.strike, args.T, t));
            std::unique_ptr<PDESolver> euroSol;
            if (solverNm == "explicit") euroSol = std::make_unique<ExplicitSolver>(*euro, market, args.M, args.N, args.Smax);
            else if (solverNm == "implicit") euroSol = std::make_unique<ImplicitSolver>(*euro, market, args.M, args.N, args.Smax);
            else euroSol = std::make_unique<CrankNicolsonSolver>(*euro, market, args.M, args.N, args.Smax);

            // Run both to fill grids
            euroSol->price();
            outSol->price();

            // Difference grid: in = euro - out (store in outSol grid for downstream)
            const Grid& G_e = euroSol->getGrid();
            Grid& G_o = outSol->getGrid();
            int Mx = G_o.getSizeS();
            int Nx = G_o.getSizeT();
            for (int j = 0; j <= Nx; ++j) {
                for (int i = 0; i <= Mx; ++i) {
                    double val = G_e.get(i,j) - G_o.get(i,j);
                    G_o.set(i,j,val);
                }
            }
            // Use outSol and outOpt as our solver/option carrying the in-grid
            opt = std::move(outOpt);
            solver = std::move(outSol);
        }
    } else if (optClass == "digital") {
        OptionType t = (optType == "put") ? OptionType::Put : OptionType::Call;
        opt = std::unique_ptr<Option>(new DigitalOption(args.strike, args.T, t, args.digitalCash));
        if (solverNm == "explicit") solver = std::make_unique<ExplicitSolver>(*opt, market, args.M, args.N, args.Smax);
        else if (solverNm == "implicit") solver = std::make_unique<ImplicitSolver>(*opt, market, args.M, args.N, args.Smax);
        else solver = std::make_unique<CrankNicolsonSolver>(*opt, market, args.M, args.N, args.Smax);
    } else {
        std::cerr << "Unknown option class: " << args.optionClass << "\n";
        print_usage();
        return 1;
    }

    double px = 0.0;
    auto t_start = std::chrono::high_resolution_clock::now();
    
    // For knock-in options, grid already contains the parity value (V_in = V_euro - V_out)
    // Do NOT recalculate with solver->price() as that would overwrite the modified grid
    bool isKnockIn = (optClass == "barrier" && 
                      (args.barrierType == "downin" || args.barrierType == "upin"));
    
    if (isKnockIn) {
        // Knock-in grid is already computed via parity: V_in = V_euro - V_out
        // Just extract the price at strike from the pre-computed grid
        const Grid& grid = solver->getGrid();
        double dS = grid.getdS();
        int iStrike = static_cast<int>(opt->getStrike() / dS);
        px = grid.get(iStrike, 0);
    } else {
        // For all other options, solve normally
        px = solver->price();
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    std::cout << "Price at strike: " << px << "\n";
    std::cout << "Time: " << std::fixed << std::setprecision(2) << elapsed_ms << " ms\n";

    // Export grid
    write_grid_csv(solver->getGrid(), opt->getMaturity(), args.out);
    std::cout << "Wrote grid CSV: " << args.out << "\n";
    // Early exercise boundary
    std::string boundaryPath;
    if (!args.earlyBoundaryOut.empty()) {
        if (write_early_boundary_csv(solver->getGrid(), *opt, opt->getMaturity(), args.earlyBoundaryOut)) {
            boundaryPath = args.earlyBoundaryOut;
            std::cout << "Wrote early boundary: " << boundaryPath << "\n";
        }
    }

    // Plot directly
    if (args.plot || args.viewer) {
        try {
            auto up = [](std::string s){ for (auto &c: s) c = std::toupper(c); return s; };
            std::string ttl = "Option 3D graph — " + args.optionClass + " " + args.optionType + " — " + up(args.solver)
                              + " | K=" + std::to_string(args.strike)
                              + ", T=" + std::to_string(args.T)
                              + ", r=" + std::to_string(args.r)
                              + ", σ=" + std::to_string(args.sigma);
            double clampK = args.clampK;
            if (clampK <= 0 && optClass == "european" && (optType == "call" || optType == "put")) {
                clampK = args.strike;
            }
            std::string overlay = (args.overlayEarlyBoundary ? boundaryPath : "");
            // Determine output mode
            std::string fmt = args.viewer ? std::string("") : args.plotFormat;
            // Default output filename if saving
            std::string outImg = args.plotOut;
            if (outImg.empty() && !fmt.empty()) {
                outImg = (fmt == "pdf") ? std::string("option_surface.pdf") : std::string("option_surface.png");
            }
            plot_surface_gnuplot(
                solver->getGrid(), opt->getMaturity(), ttl,
                args.spin, clampK, overlay,
                200, 200,
                fmt, outImg,
                args.plotWidth, args.plotHeight,
                args.plotDpi
            );
        } catch (const std::exception& e) {
            std::cerr << "Plotting failed: " << e.what() << "\n";
        }
    }

    // Optional: print space-delimited plot data to stdout (first N lines)
    if (args.printPlotData) {
        print_plot_data_stdout(solver->getGrid(), std::max(1, args.printCount));
    }
    return 0;
}



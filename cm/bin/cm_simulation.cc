#include "CLog.hh"
#include "CMSimulation.hh"
#include "CmdLineConfig.hh"
#include "Sources/MultiPointSource.hh"
#include "Sources/PlanarSource.hh"
#include "Sources/PointSource.hh"
#include <TSystem.h>
#include <iostream>

int main(int argc, char** argv) {

  CmdLineOption opt_output("Output", "-o",
                           "Output file (string), default: simulation.root",
                           "simulation.root");
  CmdLineOption opt_iter("Iterations", "-n",
                         "Number of iterations, default: 20 (integer)", 20);
  CmdLineOption opt_source_p("Point", "-point",
                             "Point-type source, default: YES");
  CmdLineOption opt_source_mp("MultiPoint", "-multipoint",
                              "Multi Point-type source, default: NO");
  CmdLineOption opt_source_pl("Planar", "-planar",
                              "Planar-type source, default: NO");
  CmdLineOption opt_plane_dim(
      "Plane", "-plane",
      "Detector plane, 6 values required, default: 1:0:0:600:300:300", 0, 0);
  CmdLineOption opt_mask_dim(
      "Mask", "-mask",
      "Mask plane, 6 values required, default: 1:0:0:500:300:300", 0, 0);

  CmdLineConfig::instance()->ReadCmdLine(argc, argv);

  PositionalArgs args = CmdLineOption::GetPositionalArguments();
  if (args.size() < 2) {
    spdlog::error("{} requires 2 arguments: mask source", argv[0]);
    abort();
  }

  spdlog::set_level(spdlog::level::info);
  TString path =
      TString(gSystem->Getenv("CC6DIR")) + "/share/ComptonCamera6/masks/";
  TString fullname = path + "hMURA" + args[0] + ".root";
  TFile* maskfile = new TFile(fullname, "READ");
  if (maskfile == nullptr) {
    spdlog::error("File with mask: {} does not exist, exit...",
                  fullname.Data());
    abort();
  }
  TH2F* h = TString(maskfile->GetName()).Contains("2d")
                ? (TH2F*)maskfile->Get("hMURA2d")
                : (TH2F*)maskfile->Get("hMURA1d");

  Float_t da = 1., db = 0., dc = 0., dd = 600., de = 300., df = 300.;
  Float_t ma = 1., mb = 0., mc = 0., md = 500., me = 300., mf = 300.;

  if (opt_plane_dim.GetArraySize() == 6) {
    da = opt_plane_dim.GetDoubleArrayValue(1);
    db = opt_plane_dim.GetDoubleArrayValue(2);
    dc = opt_plane_dim.GetDoubleArrayValue(3);
    dd = opt_plane_dim.GetDoubleArrayValue(4);
    de = opt_plane_dim.GetDoubleArrayValue(5);
    df = opt_plane_dim.GetDoubleArrayValue(6);
  } else if (opt_plane_dim.GetArraySize() != 0) {
    spdlog::error("Detector plane - 6 parameters required, {} given",
                  opt_plane_dim.GetArraySize());
    abort();
  }
  if (opt_mask_dim.GetArraySize() == 6) {
    ma = opt_mask_dim.GetDoubleArrayValue(1);
    mb = opt_mask_dim.GetDoubleArrayValue(2);
    mc = opt_mask_dim.GetDoubleArrayValue(3);
    md = opt_mask_dim.GetDoubleArrayValue(4);
    me = opt_mask_dim.GetDoubleArrayValue(5);
    mf = opt_mask_dim.GetDoubleArrayValue(6);
  } else if (opt_mask_dim.GetArraySize() != 0) {
    spdlog::error("Mask plane - 6 parameters required, {} given",
                  opt_mask_dim.GetArraySize());
    abort();
  }
  printf("Detector plane: %g %g %g %g %g %g\n", da, db, dc, dd, de, df);
  printf("Mask plane    : %g %g %g %g %g %g\n", ma, mb, mc, md, me, mf);
  printf("Mask type     : %s\n", args[0].Data());
  printf("Source file   : %s\n", args[1].Data());
  printf("Source type   : [%c] point  [%c] multi-point  [%c] planar\n",
         !opt_source_pl.GetFlagValue() and !opt_source_mp.GetFlagValue() ? 'x'
                                                                         : ' ',
         !opt_source_pl.GetFlagValue() and opt_source_mp.GetFlagValue() ? 'x'
                                                                        : ' ',
         opt_source_pl.GetFlagValue() ? 'x' : ' ');

  Source* src = nullptr;
  if (opt_source_pl.GetFlagValue())
    src = new PlanarSource(args[1]);
  else if (opt_source_mp.GetFlagValue())
    src = new MultiPointSource(args[1]);
  else if (opt_source_p.GetFlagValue())
    src = new PointSource(args[1]);
  else
    src = new PointSource(args[1]);
  src->Print();

  DetPlane detector(da, db, dc, dd, de, df, "detector");
  Mask mask(ma, mb, mc, md, me, mf, h, "mask");

  CMSimulation sim(src, &mask, &detector);
  sim.RunSimulation(CmdLineOption::GetIntValue("Iterations"));
  sim.Write(opt_output.GetStringValue());

  return 1;
}

#include "run.h"

void instance(Problem &problem, int argc, char *argv[])
{
  string instance(argv[1]);

  if (instance == "SIZES")
  {
    cout << "SIZES" << argv[2] << endl;
    problem.sizes(stoi(argv[2]));
    problem.enforce_ccr(1e4);
  }
  if (instance == "DCAP")
  {
    cout << "DCAP_" << argv[2] << '_' << argv[3] << '_' << argv[4] << '_' << argv[5] << ' ' << argv[6] << '\n';
    problem.dcap(stoi(argv[2]), stoi(argv[3]), stoi(argv[4]), stoi(argv[5]), stoi(argv[6]));
  }
  if (instance == "SSV")
  {
    cout << "SSV_" << argv[2] << '_' << argv[3] << '_' << argv[4] << '_' << argv[5] << '\n';
    problem.ssv95(stoi(argv[2]), stoi(argv[3]), stoi(argv[4]), stoi(argv[5]));
    //problem.enforce_ccr(1e4);
  }
  if (instance == "CAROE")
  {
    cout << "CAROE INSTANCE, S = " << argv[2] << '\n';
    problem.caroe( stoi(argv[2]));
    problem.enforce_ccr(1e4);
  }
  if (instance == "LARGE")
  {
    cout << "LARGE SSV, n1 = " << argv[2] << " n2 = " << argv[3] << " m2 = " << argv[4] <<
            " S = " << argv[5] << " fs_cont: " << argv[6] << " ss_bin: " << argv[7] << '\n';

    problem.ssv_large(stoi(argv[2]),
                      stoi(argv[3]),
                      stoi(argv[4]),
                      stoi(argv[5]),
                      stoi(argv[6]),
                      stoi(argv[7]));

  }
}


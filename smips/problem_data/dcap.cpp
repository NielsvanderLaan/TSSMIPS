#include "problem.h"

void Problem::dcap(size_t nResources, size_t nClients, size_t nPeriods, size_t S)
{
  d_fix_rec = false;
  d_L = 0;

  size_t nCaps = nResources * nPeriods;
  size_t combs= nClients * nResources * nPeriods;

  d_m1 = nCaps;
  d_n1 = 2 * nCaps;
  d_p1 = nCaps;
  d_m2 = nCaps + nClients * nPeriods;
  d_n2 = nResources * nClients * nPeriods + nClients * nPeriods;
  d_p2 = d_n2;
  d_S = S;

  d_fs_leq = d_m1; d_fs_geq = 0;
  d_ss_leq = nCaps; d_ss_geq = 0;

  d_l1 = vector<double> (d_n1, 0.0);
  d_u1 = vector<double> (d_n1, 1.0);
  d_l2 = vector<double> (d_n2, 0.0);
  d_u2 = vector<double> (d_n2, 1.0);


  d_Amat = vector<vector<double>> ();
  for (size_t con = 0; con != d_m1; ++con)
  {
    vector<double> row (d_n1, 0.0);
    row[con] = -1.0;
    row[d_p1 + con] = 1.0;
    d_Amat.push_back(row);
  }
  d_b = vector<double> (d_m1, 0.0);



  d_Tmat = vector<vector<double>> ();
  for (size_t res = 0; res != nResources; ++res)
  {
    vector<double> row (d_n1);
    for (size_t period = 0; period != nPeriods; ++period)
    {
      row[d_p1 + res * nPeriods + period] = -1.0;
      d_Tmat.push_back(row);
    }
  }
  for (size_t con = 0; con != nClients * nPeriods; ++con)
    d_Tmat.push_back(vector<double>(d_n1,0.0));

  vector<double> rhs (d_m2, 1.0);
  for (auto con = 0; con != nCaps; ++con)
    rhs[con] = 0.0;
  d_omega = vector<vector<double>> (S, rhs);
  d_q = vector<double> (d_n2);
  d_Wmat = vector<vector<double>> (nCaps, vector<double>(d_n2));

  for (size_t client = 0; client != nClients ; ++client)
  {
    for (size_t period = 0; period != nPeriods; ++period)
    {
      vector<double> row (d_n2);
      row[combs + client * nPeriods + period] = 1.0;
      for (size_t resource = 0; resource != nResources; ++resource)
        row[resource * nClients * nPeriods  + client * nPeriods + period] = 1.0;
      d_Wmat.push_back(row);
    }
  }


  ifstream lpFile;
  string id = to_string(nResources)+ to_string(nClients) + to_string(nPeriods);
  string base = "dcap/" + id + "/dcap" + id + "_" + to_string(S);
  string filename = base + ".lp";

  lpFile.open(filename);
  string line;
  string previous;

  d_c = vector<double> (d_n1);

  while (lpFile >> line)
  {
    if (line == "Subject")
     break;

    if (line[0] == 'u')
    {
      size_t resource = stoi(line.substr(2,1)) - 1;
      size_t period = stoi(line.substr(4, 1)) - 1;
      d_c[resource * nPeriods + period] = stod(previous);
    }

    if (line[0] == 'x')
    {
      size_t resource = stoi(line.substr(2,1)) - 1;
      size_t period = stoi(line.substr(4, 1)) - 1;
      d_c[nCaps + resource * nPeriods + period ] = stod(previous);
    }

    if (line[0] == 'y')
    {
      size_t resource = stoi(line.substr(2,1)) - 1;
      size_t client = stoi(line.substr(4, 1)) - 1;
      size_t period = stoi(line.substr(6, 1)) - 1;
      d_q[resource * nClients * nPeriods  + client * nPeriods + period] = stod(previous);
    }

    if (line[0] == 'z')
    {
      size_t client = stoi(line.substr(2, 1)) - 1;
      size_t period = stoi(line.substr(4, 1)) - 1;
      d_q[combs + client * nPeriods + period] = stod(previous);
    }

    previous = line;
  }
  lpFile.close();

  d_q_omega = vector<vector<double>> (S, d_q);
  d_probs = vector<double> (S,1.0 / S);

  ifstream stoFile;
  filename = base + ".sto";
  stoFile.open(filename);
  d_W_omega = vector<vector<vector<double>>> (S, d_Wmat);
  size_t s = -1;
  size_t row, col;
  bool value_ahead = false;

  while (stoFile >> line)
  {
    if (value_ahead)
    {
      d_W_omega[s][row][col] = stod(line);
      value_ahead = false;
    }

    if (line == "SC")
      ++s;

    if (line[0] == 'y')
    {
      size_t resource = stoi(line.substr(2,1)) - 1;
      size_t client = stoi(line.substr(4, 1)) - 1;
      size_t period = stoi(line.substr(6, 1)) - 1;
      col = resource * nClients * nPeriods  + client * nPeriods + period;
    }

    if (line.substr(0, 3) == "dem")
    {
      size_t resource = stoi(line.substr(4,1)) - 1;
      size_t period = stoi(line.substr(6, 1)) - 1;
      row = resource * nPeriods + period;
      value_ahead = true;
    }

  }
  stoFile.close();


  //d_p1 = d_n1;

}
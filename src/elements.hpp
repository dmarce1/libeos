/*
 * elements.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: dmarce1
 */


#include "real.hpp"
#include "physcon.hpp"

struct element_t {
	int Z;
	real A;
	std::vector<int> g;
	std::vector<real> e_i;
	element_t(int z, real a, std::vector<int>&& g_, std::vector<real>&& e_i_) :
			Z(z), A(a), g(std::move(g_)), e_i(std::move(e_i_)) {
		for (int i = 0; i <= Z; i++) {
			e_i[i] *= eV;
		}
	}
	real saha(int i, real T) const {
		return real(g[i + 1]) * std::exp(-e_i[i + 1] / (kb * T)) / real(g[i]);
	}
};

element_t elements[28] = {
{1,  1.0079,{{2, 1}},                                                             {{0.0, 13.59844}}},
{2,  4.0026,{{1, 2, 1}},                                                          {{0.0, 24.58741, 54.41778}}},
{3,  6.9410,{{2, 1, 2,  1}},                                                      {{0.0,  5.39172, 75.64018, 122.45429}}},
{4,  9.0122,{{1, 2, 1,  2, 1}},                                                   {{0.0,  9.3227,  18.21116, 153.89661, 217.71865}}},
{5, 10.8110,{{6, 1, 2,  1, 2, 1}},                                                {{0.0,  8.29803, 25.15484,  37.93064, 259.37521, 340.22580}}},
{6, 12.0107,{{9, 6, 1,  2, 1, 2, 1}},                                             {{0.0, 11.26030, 24.38332,  47.8878,   64.4939,  392.087,  489.99334}}},
{7, 14.0067,{{4, 9, 6,  1, 2, 1, 2, 1}},                                          {{0.0, 14.53414, 29.6013,   47.44924,  77.4735,   97.8902, 552.0718, 667.046}}},
{8, 15.9994,{{9, 4, 9,  6, 1, 2, 1, 2,1}},                                        {{0.0, 13.61806, 35.11730,  54.9355,   77.41353, 113.8990, 138.1197, 739.29,   871.4101}}},
{9, 18.9984,{{6, 9, 4,  9, 6, 1, 2, 1,2,1}},                                      {{0.0, 17.42282, 34.97082,  62.7084,   87.1398,  114.2428, 157.1651, 185.186,  953.9112, 1103.1176}}},
{10,20.1797,{{1, 6, 9,  4, 9, 6, 1, 2,1,2,1}},                                    {{0.0, 21.5646,  40.96328,  63.45,     97.12,    126.21,   157.93,   207.2759, 239.0989, 1195.8286, 1362.1995}}},
{11,22.9897,{{2, 1, 6,  9, 4, 9, 6, 1,2,1,2,1}},                                  {{0.0,  5.13908, 47.2864,   71.6200,   98.91,    138.40,   172.18,   208.50,   264.25,    299.864,  1465.121, 1648.702}}},
{12,24.3050,{{1, 2, 1,  6, 9, 4, 9, 6,1,2,1,2,1}},                                {{0.0,  7.64624, 15.03528,  80.1437,  109.2655,  141.27,   186.76,   225.02,   265.96,    328.06,    367.50,  1761.805, 1962.6650}}},
{13,26.9815,{{6, 1, 2,  1, 6, 9, 4, 9,6,1,2,1,2,1}},                              {{0.0,  5.98577, 18.82856,  28.44765, 119.992,   153.825,  190.49,   241.76,   284.66,    330.13,    398.75,   442.00,  2085.98, 2304.1410}}},
{14,28.0855,{{9, 6, 1,  2, 1, 6, 9, 4,9,6,1,2,1,2,1}},                            {{0.0,  8.15169, 16.34585,  33.49302,  45.14181, 166.767,  205.27,   246.5,    303.54,    351.12,    401.37,   476.36,   523.42, 2437.63, 2673.182}}},
{15,30.9738,{{4, 9, 6,  1, 2, 1, 6, 9,4,9,6,1,2,1,2,1}},                          {{0.0, 10.48669, 19.7694,   30.2027,   51.4439,   65.0251, 220.421,  263.57,   309.60,    372.13,    424.4,    479.46,   560.8,   611.74, 2816.91,  3069.842}}},
{16,32.0650,{{9, 4, 9,  6, 1, 2, 1, 6,9,4,9,6,1,2,1,2,1}},                        {{0.0, 10.36001, 23.3379,   34.79,     47.222,    72.5945,  88.0530, 280.948,  328.75,    379.55,    447.5,    504.8,    564.44,  652.2,   707.01,  3223.78,  3494.1892}}},
{17,35.4530,{{6, 9, 4,  9, 6, 1, 2, 1,6,9,4,9,6,1,2,1,2,1}},                      {{0.0, 12.96764, 23.814,    39.61,     53.4652,   67.8,     97.03,   114.1958, 348.28,    400.06,    455.63,   529.28,   591.99,  656.71,  749.76,   809.40,  3658.521, 3946.2960}}},
{18,39.0983,{{1, 6, 9,  4, 9, 6, 1, 2,1,6,9,4,9,6,1,2,1,2,1}},                    {{0.0, 15.75962, 27.62967,  40.74,     59.81,     75.02,    91.009,  124.323,  143.460,   422.45,    478.69,   538.96,   618.26,  686.10,  755.74,   854.77,   918.03,  4120.8857, 4426.2296}}},
{19,39.9480,{{2, 1, 6,  9, 4, 9, 6, 1,2,1,6,9,4,9,6,1,2,1,2,1}},                  {{0.0,  4.34066, 31.63,     45.806,    60.91,     82.66,    99.4,    117.56,   154.88,    175.8174,  503.8,    564.7,    629.4,   714.6,   786.6,    861.1,    968.0,   1033.4,    4610.8,  4934.046}}},
{20,40.0780,{{1, 2, 1,  6, 9, 4, 9, 6,1,2,1,6,9,4,9,6,1,2,1,2,1}},                {{0.0,  6.11316, 11.87172,  50.9131,   67.27,     84.50,   108.78,   127.2,    147.24,    188.54,    211.275,  591.9,    657.2,   726.6,   817.6,    894.5,    974.0,   1087.0,    1157.8,  5128.8,  5469.864}}},
{21,44.9559,{{10,15,10, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}},              {{0.0,  6.5615,  12.79967,  24.75666,  73.4894,   91.65,   110.68,   138.0,    158.1,     180.03,    225.18,   249.798,  687.36,  756.7,   830.8,    927.5,   1009.0,   1094.0,    1213.0,  1287.97, 5674.8, 6033.712}}},
{22,47.8670,{{21,28,21,10, 1, 6, 9, 4,9,2,1,2,1,2,1,2,1,2,1,2,1,2,1}},            {{0.0,  6.8281,  13.5755,   27.4917,   43.2672,   99.30,   119.53,   140.8,    170.4,     192.1,     215.92,   265.07,   291.500, 787.84,  863.1,    941.9,   1044.0,   1131.0,    1221.0,  1346.0,  1425.4, 6249.0, 6625.82}}},
{23,50.9415,{{28,25,28, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}},          {{0.0,  6.7463,  14.66,     29.311,    46.709,    65.2817, 128.13,   150.6,    173.4,     205.8,     230.5,    255.7,    308.1,   336.277, 896.0,    976.0,   1060.0,   1168.0,    1260.0,  1355.0,  1486.0, 1569.6, 6851.3, 7246.12}}},
{24,51.9961,{{ 7, 6,25,28,21,10, 1, 6,9,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1}},        {{0.0,  6.7665,  16.4857,   30.96,     49.16,     69.46,    90.6349, 160.18,   184.7,     209.3,     244.4,    270.8,    298.0,   354.8,   384.168, 1010.6,   1097.0,   1185.0,    1299.0,  1396.0,  1496.0, 1634.0, 1721.4, 7481.7, 7894.81}}},
{25,54.9380,{{ 6, 7, 6, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}},      {{0.0,  7.43402, 15.63999,  33.668,    51.2,      72.4,     95.6,    119.203,  194.5,     221.8,     248.3,    286.0,    314.4,   343.6,   403.0,    435.163, 1134.7,   1224.0,    1317.0,  1437.0,  1539.0, 1644.0, 1788.0, 1879.9, 8140.6, 8571.94}}},
{26,55.8450,{{25,30,25, 6,25,28,21,10,1,6,9,2,9,2,1,2,1,2,1,2,1,2,1,2,1,2,1}},    {{0.0,  7.9024,  16.1878,   30.652,    54.8,      75.0,     99.1,    124.98,   151.06,    233.6,     262.1,    290.2,    330.8,   361.0,   392.2,    457.0,    489.256, 1266.0,    1358.0,  1456.0,  1582.0, 1689.0, 1799.0, 1950.0, 2023.0, 8828.0, 9277.69}}},
{27,58.9332,{{ 1, 1, 1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}},  {{0.0,  7.8810,  17.083,    33.50,     51.3,      79.5,    102.0,    128.9,    157.8,     186.13,    275.4,    305.0,    336.0,   379.0,   411.0,    444.0,    511.96,   546.58,   1397.2,  1504.6,  1603.0, 1735.0, 1846.0, 1962.0, 2119.0, 2219.0, 9544.1, 10012.12}}},
{28,58.6934,{{21,10,21,10,10,10,25,28,6,6,6,6,9,6,9,6,1,2,1,2,1,2,1,2,1,2,1,2,1}},{{0.0,  7.6398,  18.16884,  35.19,     54.9,      76.06,   108.0,    133.0,    162.0,     193.0,     224.6,    321.0,    352.0,   384.0,   430.0,    464.0,    499.0,    571.08,    607.06, 1541.0,  1648.0, 1756.0, 1894.0, 2011.0, 2131.0, 2295.0, 2399.2, 10288.8, 10775.40}}}
};
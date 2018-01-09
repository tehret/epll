/*
 * Original work: Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EPLL_H_INCLUDED
#define EPLL_H_INCLUDED

#define NTHREAD 4
#define PI 3.1415926535

#include <vector>
#include "../Utilities/LibVideoT.hpp"

struct Model{
	float logweight;
	std::vector<float> eigVects;
	std::vector<float> eigVals;
	std::vector<float> invSqrtCov;
	float logdet;
	int r;
};

void EPLLhalfQuadraticSplit(Video<float>& noiseI, Video<float>& finalI, Video<float>& origI, bool partialPSNR, float noiseSD, int patchsize, int patchsizeChannels, std::vector<float> betas, int T, int pas, std::vector<Model>& models);

void aprxMAPGMM(Video<float>& noiseI, Video<float>& tempI, float sigma, int ps, int psc, int step, std::vector<Model>& models);

void loggausspdf(std::vector<float>& patches, int dim, int nbP, Model& model, float csta, std::vector<float>& result);

#endif

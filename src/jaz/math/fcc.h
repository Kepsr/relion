#ifndef RELION_TOMO_FCC_H
#define RELION_TOMO_FCC_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/particle_set.h>


class FCC
{
	public:
		
		static BufferedImage<double> compute(
				ParticleSet* dataSet,
				const std::vector<int>& partIndices,
				const Tomogram& tomogram,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				bool flip_value,
				int num_threads);
		
		static BufferedImage<double> compute3(
				ParticleSet* dataSet,
				const std::vector<int>& partIndices,
				const Tomogram& tomogram,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				bool flip_value,
				int num_threads);
		
		static BufferedImage<double> divide(
				const BufferedImage<double>& fcc3);
};

#endif

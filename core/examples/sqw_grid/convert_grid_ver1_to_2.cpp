/**
 * converts version 1 grid data to version 2
 * @author Tobias Weber <tweber@ill.fr>
 * @date 05-jan-20
 * @license GPLv2
 */

/*
 * Grid version 2 format
 *
 * Header format:
 *     8 bytes (std::size_t): offset of index block
 *     3*8 bytes (double): data dimensions: hmin, hmax, hstep
 *     3*8 bytes (double): data dimensions: kmin, kmax, kstep
 *     3*8 bytes (double): data dimensions: lmin, lmax, lstep
 *     x bytes: metadata header
 *
 * <data block>
 * <index block>
 *
 * Data block format:
 *     repeat for each reduced wave vector q:
 *         4 bytes (unsigned int): number of dispersion branches
 *             repeat (number of branches times):
 *                 8 bytes (double): energy
 *                 8 bytes (double): dynamical structure factor
 *
 * Index block format:
 *     repeat for each (h,k,l) coordinate:
 *         8 bytes (std::size_t): offset into data block
 */

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

using t_float_src = double;
using t_float_dst = double;


int main()
{
	const t_float_dst eps = 1e-8;


	std::string filenameIdx = "grid_ver1.idx";
	std::string filenameDat = "grid_ver1.bin";

	std::string filenameNewDat = "grid_ver2.bin";


	// dimensions of version 1 grid (TODO: load from config file)
	t_float_dst hmin = -0.096;
	t_float_dst hmax = 0.096;
	t_float_dst hstep = 0.002;

	t_float_dst kmin = -0.096;
	t_float_dst kmax = 0.096;
	t_float_dst kstep = 0.002;

	t_float_dst lmin = -0.096;
	t_float_dst lmax = 0.096;
	t_float_dst lstep = 0.002;


	std::cout << "Converting data file ..." << std::endl;

	std::ifstream ifDat(filenameDat);
	std::ofstream ofDat(filenameNewDat);

	std::unordered_map<std::size_t, std::size_t> mapIndices;

	std::size_t removedBranches = 0;


	// write a dummy index file offset at the beginning (to be filled in later)
	std::size_t idx_offs = 0;
	ofDat.write((char*)&idx_offs, sizeof(idx_offs));


	// write data dimensions
	ofDat.write((char*)&hmin, sizeof(hmin));
	ofDat.write((char*)&hmax, sizeof(hmax));
	ofDat.write((char*)&hstep, sizeof(hstep));

	ofDat.write((char*)&kmin, sizeof(kmin));
	ofDat.write((char*)&kmax, sizeof(kmax));
	ofDat.write((char*)&kstep, sizeof(kstep));

	ofDat.write((char*)&lmin, sizeof(lmin));
	ofDat.write((char*)&lmax, sizeof(lmax));
	ofDat.write((char*)&lstep, sizeof(lstep));


	// header metadata
	ofDat << "takin_grid_data_ver2|title:TEST|author:tweber@ill.fr|date:5/feb/2020";


	while(1)
	{
		std::size_t idx = ifDat.tellg();
		std::size_t idxnew = ofDat.tellp();

		unsigned int numBranches = 0;
		ifDat.read((char*)&numBranches, sizeof(numBranches));
		if(ifDat.gcount() != sizeof(numBranches) || ifDat.eof())
			break;

		mapIndices.insert(std::make_pair(idx, idxnew));
		//std::cout << "index " << std::hex << idx << " -> " << idxnew << "\n";


		// write placeholder
		unsigned int numNewBranches = 0;
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));


		for(unsigned int branch=0; branch<numBranches; ++branch)
		{
			// in case more than one weigh factor is stored in original grid:
			//t_float_src vals[4] = {0, 0, 0, 0};

			t_float_src vals[2] = { 0, 0 };
			ifDat.read((char*)vals, sizeof(vals));

			t_float_dst E = vals[0];
			if(std::abs(E) < eps)
				E = t_float_dst(0);

			//t_float_dst w = std::abs(vals[1])+std::abs(vals[2])+std::abs(vals[3]);
			t_float_dst w = std::abs(vals[1]);

			if(w >= eps)
			{
				t_float_dst newvals[2] = { E, w };
				ofDat.write((char*)newvals, sizeof(newvals));

				++numNewBranches;
			}
			else
			{
				++removedBranches;
			}
		}


		// seek back and write real number of branches
		std::size_t lastIdx = ofDat.tellp();
		ofDat.seekp(idxnew, std::ios_base::beg);
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));
		ofDat.seekp(lastIdx, std::ios_base::beg);
	}

	ifDat.close();

	// update index at beginning
	idx_offs = ofDat.tellp();
	ofDat.seekp(0, std::ios_base::beg);
	ofDat.write((char*)&idx_offs, sizeof(idx_offs));
	ofDat.seekp(idx_offs, std::ios_base::beg);


	std::cout << removedBranches << " branches removed (w < eps)." << std::endl;


	std::cout << "\nConverting index file ..." << std::endl;

	std::ifstream ifIdx(filenameIdx);
	// write index block in the same file as data
	std::ofstream &ofIdx = ofDat;


	while(1)
	{
		std::size_t idx = 0;

		ifIdx.read((char*)&idx, sizeof(idx));
		if(ifIdx.gcount() != sizeof(idx) || ifIdx.eof())
			break;

		auto iter = mapIndices.find(idx);
		if(iter == mapIndices.end())
		{
			std::cerr << "Error: Index " << std::hex << idx << " was not found." << std::endl;
			continue;
		}

		std::size_t newidx = iter->second;
		ofIdx.write((char*)&newidx, sizeof(newidx));
	}


	return 0;
}

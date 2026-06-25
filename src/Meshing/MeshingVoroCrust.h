///////////////////////////////////////////////////////////////////////////////////////////////
//                                VOROCRUST-MESHING 1.0                                      //
// Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).        //
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain  //
// rights in this software.                                                                  //
//                                                                                           //
// Redistribution and use in source and binary forms, with or without modification, are      //
// permitted provided that the following conditions are met:                                 //  
//                                                                                           //
// 1. Redistributions of source code must retain the above copyright notice, this list of    //
// conditions and the following disclaimer.                                                  //
//                                                                                           //
// 2. Redistributions in binary form must reproduce the above copyright notice, this list    //
// of conditions and the following disclaimer in the // documentation and/or other materials //
// provided with the distribution.                                                           //
//                                                                                           //
// 3. Neither the name of the copyright holder nor the names of its contributors may be      //
// used to endorse or promote products derived from this software without specific prior     //
// written permission.                                                                       //
//-------------------------------------------------------------------------------------------//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY       //
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF   //
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE//
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, //
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        //
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    //
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR  //
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                              //
///////////////////////////////////////////////////////////////////////////////////////////////
//                                     Author                                                //
//                                Mohamed S. Ebeida                                          //
//                                msebeid@sandia.gov                                         //
///////////////////////////////////////////////////////////////////////////////////////////////
//  MeshingVoroCrust.h                                            Last modified (08/03/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _VOROCRUST_H_
#define _VOROCRUST_H_

#include <string>

#include <cmath>
#include <list>
#include <vector>
#include <utility>
#include <stack>
#include <map>
#include <memory>
#include <set>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <iomanip>

#include <time.h>

#include "Config.h"
#include "MeshingConsole.h"
#include "MeshingOptionParser.h"
#include "MeshingTimer.h"

#if defined USE_OPEN_MP
#include <omp.h>
#endif


#ifndef DBL_MAX
#define DBL_MAX  1.7976931348623158e+308 // max double value
#endif

#ifndef PI
#define PI  3.141592653589793
#endif

// extern std::shared_ptr<VoroCrust_OptionParser> options;

extern MeshingConsole vcm_cout;


#include "MeshingSmartTree.h"
#include "MeshingRandomSampler.h"

#include "MeshingVoroCrustSampler.h"
#include "MeshingVoronoiMesher.h"
#include "MeshingVoroCrustObserver.h"

struct MeshingVoroCrustOptions
{
	std::vector<std::string> input_mesh_files;   ///< 输入网格文件列表（OBJ 格式），多个文件会自动合并
	std::string sizing_file;                     ///< 可选的尺寸场文件（CSV 格式），空字符串表示不使用
	std::string seeds_file;                      ///< 种子点文件（CSV 格式），仅在 generate_mesh_from_seeds 模式下使用
	std::string input_filename;                  ///< 可选的 .in 输入文件名；若设置，则优先解析该文件

	double r_min = 0.0;                          ///< 球体最小半径，0 表示根据模型自动处理干净输入
	double r_max = DBL_MAX;                      ///< 球体最大半径
	double Lip_const = 0.1;                      ///< Lipschitz 常数，控制尺寸场变化速率
	double vc_ang_tol = 20.0;                    ///< 特征检测/平滑的角度阈值（度）
	int num_threads = 0;                         ///< OpenMP 线程数；0 表示启用 OpenMP 时自动检测，否则使用 1

	bool generate_vcg_file = false;              ///< 是否生成 VCG 网格文件
	bool generate_exodus_file = false;           ///< 是否生成 Exodus 网格文件
	bool generate_monitoring_points = false;     ///< 是否生成监测点文件 monitoring_points.csv
	bool impose_monitoring_points = false;       ///< 是否从 monitoring_points.csv 强制插入监测点
	bool generate_mesh_from_seeds = false;       ///< 是否跳过采样，直接从 seeds_file 生成网格

	double feature_TOL = 0.0;                    ///< 特征容差，保留给未来特征处理使用
	size_t num_loop_ref = 0;                     ///< 循环细分/细化次数，保留给未来使用
	double ref_ang_tol = 20.0;                   ///< 细化角度阈值（度），保留给未来使用
};


class MeshingVoroCrust
{

public:

	MeshingVoroCrust(std::string input_filename);
	MeshingVoroCrust(const MeshingVoroCrustOptions& options);

	~MeshingVoroCrust();

	int execute();

	const MeshingVoroCrustOptions& get_options() const { return options; }

	void set_observer(MeshingVoroCrustObserver* observer) { _observer = observer; }

private:

	int parse_input_file(MeshingVoroCrustOptions& opts);
	int execute_with_options(const MeshingVoroCrustOptions& opts);

	int execute_vc(int num_threads, size_t num_points, double** points, size_t num_faces, size_t** faces, MeshingSmartTree* input_sizing_function,
		                double rmin, double rmax, double Lip, double smooth_angle_threshold, bool generate_vcg_file,bool generate_exodus_file,
		                bool generate_monitoring_points, bool impose_monitoring_points);

	int generate_interior_seeds(MeshingSmartTree* seeds_tree, MeshingSmartTree* surface_spheres_tree, MeshingSmartTree* edge_spheres_tree, MeshingSmartTree* corner_spheres_tree,
		                        MeshingSmartTree* sz_function_tree, int num_threads, double Lip, double rmax,
		                        bool impose_monitoring_points, bool generate_monitoring_points,
		                        size_t& num_seeds, double*& seeds, size_t*& seeds_region_id, double*& seeds_sizing);

	int generate_explicit_mesh(int num_threads, size_t num_seeds, double* seeds, size_t* seeds_region_id, double* seeds_sizing,
		                       size_t num_surface_seeds, bool generate_vcg_file, bool generate_exodus_file);

	int generate_mesh_from_seeds_file(std::string file_name, int num_threads, bool generate_vcg_file, bool generate_exodus_file);

	int read_seeds_csv(std::string file_name, size_t& num_seeds, double*& seeds, size_t*& seeds_region_id, double*& seeds_sizing);

	int get_tokens(std::string line, char separator, std::vector<std::string>& tokens);


	MeshingRandomSampler _rsampler;
	MeshingMemoryHandler _memo;

	MeshingVoroCrustOptions options;
	MeshingVoroCrustObserver* _observer = 0;

};

#endif

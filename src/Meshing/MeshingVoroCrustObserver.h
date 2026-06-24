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
// of conditions and the following disclaimer in the documentation and/or other materials     //
// provided with the distribution.                                                           //
//                                                                                           //
// 3. Neither the name of the copyright holder nor the names of its contributors may be      //
// used to endorse or promote products derived from this software without specific prior     //
// written permission.                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _VOROCRUST_OBSERVER_H_
#define _VOROCRUST_OBSERVER_H_

#include "MeshingSmartTree.h"

#include <cstddef>
#include <string>
#include <vector>

struct MeshingVoroCrustTreePointSnapshot
{
	std::vector<double> coords;
	std::vector<double> normal;
	std::vector<size_t> attrib;
};

struct MeshingVoroCrustTreeSnapshot
{
	std::string name;
	size_t num_dim = 0;
	std::vector<MeshingVoroCrustTreePointSnapshot> points;
};

struct MeshingVoroCrustFeatureSnapshot
{
	MeshingVoroCrustTreeSnapshot surface_point_cloud;
	MeshingVoroCrustTreeSnapshot edge_point_cloud;
	std::vector<size_t> sharp_corners;
	std::vector<std::vector<size_t> > sharp_edges;
};

struct MeshingVoroCrustSphereSnapshot
{
	MeshingVoroCrustTreeSnapshot corner_spheres;
	MeshingVoroCrustTreeSnapshot edge_spheres;
	MeshingVoroCrustTreeSnapshot surface_spheres;
	MeshingVoroCrustTreeSnapshot surf_ext_seeds;
	MeshingVoroCrustTreeSnapshot surf_int_seeds;
};

struct MeshingVoroCrustSurfaceSeedSnapshot
{
	MeshingVoroCrustTreeSnapshot interface_spheres;
	MeshingVoroCrustTreeSnapshot seeds;
	size_t num_subregions = 0;
};

struct MeshingVoroCrustFlatSeedSnapshot
{
	std::string name;
	size_t num_seeds = 0;
	std::vector<double> xyz;
	std::vector<double> sizing;
	std::vector<size_t> region_id;
};

struct MeshingVoroCrustMeshSnapshot
{
	std::string name;
	size_t num_vertices = 0;
	std::vector<double> vertices;
	std::vector<std::vector<size_t> > faces;
	std::vector<size_t> seeds_region_id;
	std::vector<size_t> cell_num_faces;
	std::vector<std::vector<size_t> > cell_faces;
};

class MeshingVoroCrustObserver
{
public:
	virtual ~MeshingVoroCrustObserver() = default;

	virtual void on_features(const MeshingVoroCrustFeatureSnapshot&) {}
	virtual void on_spheres(const MeshingVoroCrustSphereSnapshot&) {}
	virtual void on_surface_seeds(const MeshingVoroCrustSurfaceSeedSnapshot&) {}
	virtual void on_seed_array(const MeshingVoroCrustFlatSeedSnapshot&) {}
	virtual void on_mesh(const MeshingVoroCrustMeshSnapshot&) {}
};

inline MeshingVoroCrustTreeSnapshot meshing_vorocrust_snapshot_tree(const std::string& name, MeshingSmartTree* tree, size_t num_dim)
{
	MeshingVoroCrustTreeSnapshot snapshot;
	snapshot.name = name;
	if (tree == 0) return snapshot;
	if (num_dim == 0) num_dim = tree->get_num_point_dimensions();
	snapshot.num_dim = num_dim;

	size_t num_points = tree->get_num_tree_points();
	snapshot.points.reserve(num_points);
	for (size_t ipoint = 0; ipoint < num_points; ++ipoint)
	{
		MeshingVoroCrustTreePointSnapshot point;
		double* x = tree->get_tree_point(ipoint);
		double* normal = tree->get_tree_point_normal(ipoint);
		size_t* attrib = tree->get_tree_point_attrib(ipoint);

		if (x != 0) point.coords.assign(x, x + num_dim);
		if (normal != 0) point.normal.assign(normal, normal + num_dim);
		if (attrib != 0) point.attrib.assign(attrib, attrib + attrib[0]);
		snapshot.points.push_back(point);
	}
	return snapshot;
}

inline MeshingVoroCrustFlatSeedSnapshot meshing_vorocrust_snapshot_seed_array(
	const std::string& name, size_t num_seeds, double* seeds, size_t* seeds_region_id, double* seeds_sizing)
{
	MeshingVoroCrustFlatSeedSnapshot snapshot;
	snapshot.name = name;
	snapshot.num_seeds = num_seeds;
	if (num_seeds == 0) return snapshot;

	if (seeds != 0) snapshot.xyz.assign(seeds, seeds + num_seeds * 3);
	if (seeds_sizing != 0) snapshot.sizing.assign(seeds_sizing, seeds_sizing + num_seeds);
	if (seeds_region_id != 0) snapshot.region_id.assign(seeds_region_id, seeds_region_id + num_seeds);
	return snapshot;
}

inline MeshingVoroCrustMeshSnapshot meshing_vorocrust_snapshot_mesh(
	const std::string& name,
	size_t num_vertices,
	double* vertices,
	size_t num_faces,
	size_t** faces,
	size_t num_seeds,
	size_t* seeds_region_id,
	size_t* cell_num_faces,
	size_t** cell_faces)
{
	MeshingVoroCrustMeshSnapshot snapshot;
	snapshot.name = name;
	snapshot.num_vertices = num_vertices;
	if (vertices != 0) snapshot.vertices.assign(vertices, vertices + num_vertices * 3);

	snapshot.faces.reserve(num_faces);
	for (size_t iface = 0; iface < num_faces; ++iface)
	{
		size_t* face = faces[iface];
		size_t num_entries = face[0] + 3;
		snapshot.faces.push_back(std::vector<size_t>(face, face + num_entries));
	}

	if (seeds_region_id != 0) snapshot.seeds_region_id.assign(seeds_region_id, seeds_region_id + num_seeds);
	if (cell_num_faces != 0) snapshot.cell_num_faces.assign(cell_num_faces, cell_num_faces + num_seeds);
	if (cell_faces != 0 && cell_num_faces != 0)
	{
		snapshot.cell_faces.resize(num_seeds);
		for (size_t iseed = 0; iseed < num_seeds; ++iseed)
		{
			if (cell_num_faces[iseed] == 0 || cell_faces[iseed] == 0) continue;
			snapshot.cell_faces[iseed].assign(cell_faces[iseed], cell_faces[iseed] + cell_num_faces[iseed]);
		}
	}
	return snapshot;
}

#endif

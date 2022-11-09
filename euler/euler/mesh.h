#ifndef EULER_MESH_H
#define EULER_MESH_H


#include <euler/classes.h>
#include <euler/flux.h>
#include <fstream>
#include <chrono>



// trim from left
inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }

std::string trim_all_whitespaces(std::string s) {
    std::string str = s;
    trim(str);
    std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
    str.erase(new_end, str.end()); 
    return str;
}


// Find index of element in array
int find(std::vector<std::string>& a, std::string v) {
	bool found = false;
	int n = 0;
	while ( (!found)&(n < a.size()) ) {
		if (v == a[n]) {
			found = true;
			return n;
		}
		n += 1;
	}
	return -1;
}

// Read from text file named filename
std::vector<std::string> read(const std::string filename) {
	std::string line;
	std::vector<std::string> out;
	std::ifstream f(filename);

	if (f.is_open()) {
		while (std::getline(f, line)) {
            //line.pop_back();
			out.push_back( line );
		}
        f.close();
	}

	return out;
}

// Get subvector of vector
std::vector<std::string> subvector(const std::vector<std::string>& v, const int I, const int J) {
	std::vector<std::string> out;	// Init vector
	for (int i=I; i <= J; ++i ) {
		out.push_back(v[i]);
	}
	return out;
}

// Cut string by character c
std::vector<std::string> cut_str(const std::string s, const char c) {
	std::vector<std::string> out;
	std::vector<int> I;
	// find zeros
	for (int i=0; i < s.size(); ++i) {
		if (s[i] == c) {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( s );
	} else {
		out.push_back( s.substr(0, I[0]) );
		for (int j=0; j < (I.size()-1); ++j) {
			out.push_back( s.substr(I[j]+1, I[j+1]-I[j]-1) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1) );
		}
	}
	return out;
}
std::vector<int> str_to_ints(std::string s) {
	std::vector<int> out;
	std::vector<int> I;
	// find zeros
	for (int i=0; i < s.size(); ++i) {
		if (s[i] == ' ') {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( std::stoi(s) );
	} else {
		out.push_back( std::stoi(s.substr(0, I[0])) );
		for (int j=0; j < (I.size()-1); ++j) {
			out.push_back( std::stoi(s.substr(I[j]+1, I[j+1]-I[j])) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( std::stoi(s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1)) );
		}
	}
	return out;
}
std::vector<double> str_to_floats(const std::string s) {
	std::vector<double> out;
	std::vector<int> I;
	// find zeros
	for (int i=0; i < s.size(); ++i) {
		if (s[i] == ' ') {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( std::stod(s) );
	} else {
		out.push_back( std::stod(s.substr(0, I[0])) );
		for (int j=0; j < (I.size()-1); ++j) {
			out.push_back( std::stod(s.substr(I[j]+1, I[j+1]-I[j])) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( std::stod(s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1)) );
		}
	}
	return out;
}

int find_tag(std::string tag, std::vector<std::string> s) {
    for (int i=0; i<s.size(); ++i) {
        if (s[i].find(tag) != std::string::npos) {
            return i;
        }
    }
    return s.size();
}



uint mesh::find_edge_with_nodes(const uint n0, const uint n1) {
    uint nmin, nmax;
    if (n0 < n1) {
        nmin = n0, nmax = n1;
    } else {
        nmin = n1, nmax = n0;
    }
    // Also check in bounds
    const auto tp = std::make_tuple(nmin, nmax);
    if (edgesRef.find(tp) != edgesRef.end()) {
        // Found
        return edgesRef.at(tp);
    } else {
        return -1;
    }
}

uint mesh::find_bound_with_nodes(const uint n0, const uint n1) {
    uint nmin, nmax;
    if (n0 < n1) {
        nmin = n0, nmax = n1;
    } else {
        nmin = n1, nmax = n0;
    }
    const auto tp = std::make_tuple(nmin, nmax);
    // Also check in bounds
    if (boundsRef.find(tp) != boundsRef.end()) {
        // Found
        return boundsRef.at(tp);
    } else {
        return -1;
    }
}

bool mesh::find_if_edge_in_mesh(const uint n0, const uint n1) {
    uint nmin, nmax;
    if (n0 < n1) {
        nmin = n0, nmax = n1;
    } else {
        nmin = n1, nmax = n0;
    }
    const auto tp = std::make_tuple(nmin, nmax);
    if (edgesRef.find(tp) != edgesRef.end()) {
        // Found
        return true;
    }
    // Also check in bounds
    if (boundsRef.find(tp) != boundsRef.end()) {
        // Found
        return true;
    }
    return false;
}


void mesh::convert_node_face_info() {
    // Take a mesh input from su2, and convert its info to correct fvm data
    // Currently,
    //      edges contain no face connectivity data
    //      bounds contain no face connectivity data

    // Loop over all faces
    std::vector<int> faceEdges(4);
    std::vector<int> faceBounds(4);
    for (int i=0; i<facesAreas.size(); ++i) {

		// Find edges
        faceEdges[0] = find_edge_with_nodes(facesNodes0[i], facesNodes1[i]);
        faceEdges[1] = find_edge_with_nodes(facesNodes1[i], facesNodes2[i]);
        if (faceIsTriangle[i]) {
            faceEdges[2] = find_edge_with_nodes(facesNodes2[i], facesNodes0[i]);
            faceEdges[3] = -1;
        } else {
            faceEdges[2] = find_edge_with_nodes(facesNodes2[i], facesNodes3[i]);
            faceEdges[3] = find_edge_with_nodes(facesNodes3[i], facesNodes0[i]);
        }

        // Find boundaries
        faceBounds[0] = find_bound_with_nodes(facesNodes0[i], facesNodes1[i]);
        faceBounds[1] = find_bound_with_nodes(facesNodes1[i], facesNodes2[i]);
        if (faceIsTriangle[i]) {
            faceBounds[2] = find_bound_with_nodes(facesNodes2[i], facesNodes0[i]);
            faceBounds[3] = -1;
        } else {
            faceBounds[2] = find_bound_with_nodes(facesNodes2[i], facesNodes3[i]);
            faceBounds[3] = find_bound_with_nodes(facesNodes3[i], facesNodes0[i]);
        }

        // Now update these edges with face connectivity info
        for (int ei : faceEdges) {
            if (ei != -1) {
                if (i != edgesFaces0[ei]) {
                    edgesFaces1[ei] = i;
                }
            }
        }

        // Do the same for boundaries
        for (int bi : faceBounds) {
            if (bi != -1) {
                boundsFaces[bi] = i;
            }
        }
    }
}



void mesh::compute_mesh() {

	// Compute mesh normals, areas, lengths
    // Loop over each face
	for (int i=0; i<facesAreas.size(); ++i) {
		// Evaluate face centroid position
        vec2 faceC = vec2(0., 0.);
        
        if (faceIsTriangle[i]) {
            faceC += nodes[facesNodes0[i]]/3.0;
            faceC += nodes[facesNodes1[i]]/3.0;
            faceC += nodes[facesNodes2[i]]/3.0;
        } else {
            faceC += nodes[facesNodes0[i]]/4.0;
            faceC += nodes[facesNodes1[i]]/4.0;
            faceC += nodes[facesNodes2[i]]/4.0;
            faceC += nodes[facesNodes3[i]]/4.0;
        }
        facesCenters[i]= faceC;
        facesAreas[i] = 0.;

    }

    // Evaluate edge normals and edge centers
	// Loop over each edge
	for (int i=0; i<edgesLengths.size(); ++i) {
        edgesCenters[i] = (nodes[edgesNodes1[i]] + nodes[edgesNodes0[i]])*0.5;

        auto de = nodes[edgesNodes1[i]] - nodes[edgesNodes0[i]];
        real l = norm(de);
		edgesLengths[i] = l;
		edgesNormals[i].x = -de.y/l;
        edgesNormals[i].y =  de.x/l;

        // Check normal direction
        auto dot_f_f = dot(edgesNormals[i],  edgesCenters[i] -  facesCenters[edgesFaces0[i]]);
        if (dot_f_f < 0) {
            edgesNormals[i] *= -1.;
        }
	}

    // Evaluate bound normals, bound centers
    // Loop over each bound
	for (int i=0; i<boundsLengths.size(); ++i) {
        boundsCenters[i] = (nodes[boundsNodes1[i]] + nodes[boundsNodes0[i]])*0.5;

        auto de = nodes[boundsNodes1[i]] - nodes[boundsNodes0[i]];
        real l = norm(de);
		boundsLengths[i] = l;
		boundsNormals[i].x = -de.y/l;
        boundsNormals[i].y =  de.x/l;

        // Check normal direction
        auto dot_f_f = dot(boundsNormals[i],  boundsCenters[i] -  facesCenters[boundsFaces[i]]);
        if (dot_f_f < 0) {
            boundsNormals[i] *= -1.;
        }
	}

    // Compute face area
    for (int i=0; i<edgesLengths.size(); ++i) {
        // Use divergence theorem

        facesAreas[edgesFaces0[i]] += dot( 
                edgesNormals[i],
                edgesCenters[i] - facesCenters[edgesFaces0[i]]
            )*edgesLengths[i]*0.5;
        
        facesAreas[edgesFaces1[i]] += dot( 
                edgesNormals[i]*-1.,
                edgesCenters[i] - facesCenters[edgesFaces1[i]]
            )*edgesLengths[i]*0.5;
    }
	for (int i=0; i<boundsLengths.size(); ++i) {
        // Use divergence theorem

        facesAreas[boundsFaces[i]] += dot( 
                boundsNormals[i],
                boundsCenters[i] - facesCenters[boundsFaces[i]]
            )*boundsLengths[i]*0.5;
    }

	// Loop over each face
	nodesGradWeights.resize(nodes.size());
    for (int i=0; i<nodesGradWeights.size(); ++i) {
        nodesGradWeights[i] = 0.;
    }
    for (int i=0; i<facesAreas.size(); ++i) {
        facesAreas[i] = std::abs(facesAreas[i]);

		// Also add node face gradient weights
        std::array<uint, 3> ns = {
            facesNodes0[i],
            facesNodes1[i],
            facesNodes2[i]
        };

        for (int ni : ns) {
            nodesGradWeights[ni] += 1./norm(nodes[ni] - facesCenters[i]);
		}
        if (!faceIsTriangle[i]) {
            int ni = facesNodes3[i];
            nodesGradWeights[ni] += 1./norm(nodes[ni] - facesCenters[i]);
        }
    }
}


void mesh::add_last_quad_edges() {

    int nfaces = facesNodes0.size()-1;
    if (!find_if_edge_in_mesh(facesNodes0[nfaces], facesNodes1[nfaces])) {
        edgesNormals.push_back(vec2());
        edgesCenters.push_back(vec2());
        edgesFaces0.push_back(nfaces);
        edgesFaces1.push_back(0);
        edgesNodes0.push_back(facesNodes0[nfaces]);
        edgesNodes1.push_back(facesNodes1[nfaces]);
        edgesLengths.push_back(0.);
        const uint nmin = std::min(facesNodes0[nfaces], facesNodes1[nfaces]);
        const uint nmax = std::max(facesNodes0[nfaces], facesNodes1[nfaces]);
        edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
    }
    if (!find_if_edge_in_mesh(facesNodes1[nfaces], facesNodes2[nfaces])) {
        edgesNormals.push_back(vec2());
        edgesCenters.push_back(vec2());
        edgesFaces0.push_back(nfaces);
        edgesFaces1.push_back(0);
        edgesNodes0.push_back(facesNodes1[nfaces]);
        edgesNodes1.push_back(facesNodes2[nfaces]);
        edgesLengths.push_back(0.);
        const uint nmin = std::min(facesNodes1[nfaces], facesNodes2[nfaces]);
        const uint nmax = std::max(facesNodes1[nfaces], facesNodes2[nfaces]);
        edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
    }
    if (!find_if_edge_in_mesh(facesNodes2[nfaces], facesNodes3[nfaces])) {
        edgesNormals.push_back(vec2());
        edgesCenters.push_back(vec2());
        edgesFaces0.push_back(nfaces);
        edgesFaces1.push_back(0);
        edgesNodes0.push_back(facesNodes2[nfaces]);
        edgesNodes1.push_back(facesNodes3[nfaces]);
        edgesLengths.push_back(0.);
        const uint nmin = std::min(facesNodes2[nfaces], facesNodes3[nfaces]);
        const uint nmax = std::max(facesNodes2[nfaces], facesNodes3[nfaces]);
        edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
    }
    if (!find_if_edge_in_mesh(facesNodes3[nfaces], facesNodes0[nfaces])) {
        edgesNormals.push_back(vec2());
        edgesCenters.push_back(vec2());
        edgesFaces0.push_back(nfaces);
        edgesFaces1.push_back(0);
        edgesNodes0.push_back(facesNodes3[nfaces]);
        edgesNodes1.push_back(facesNodes0[nfaces]);
        edgesLengths.push_back(0.);
        const uint nmin = std::min(facesNodes3[nfaces], facesNodes0[nfaces]);
        const uint nmax = std::max(facesNodes3[nfaces], facesNodes0[nfaces]);
        edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
    }
}


void mesh::add_last_tri_edges() {
    int nfaces = facesNodes0.size()-1;
    if (!find_if_edge_in_mesh(facesNodes0[nfaces], facesNodes1[nfaces])) {
        edgesNormals.push_back(vec2());
        edgesCenters.push_back(vec2());
        edgesFaces0.push_back(nfaces);
        edgesFaces1.push_back(0);
        edgesNodes0.push_back(facesNodes0[nfaces]);
        edgesNodes1.push_back(facesNodes1[nfaces]);
        edgesLengths.push_back(0.);
        const uint nmin = std::min(facesNodes0[nfaces], facesNodes1[nfaces]);
        const uint nmax = std::max(facesNodes0[nfaces], facesNodes1[nfaces]);
        edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
    }
    if (!find_if_edge_in_mesh(facesNodes1[nfaces], facesNodes2[nfaces])) {
        edgesNormals.push_back(vec2());
        edgesCenters.push_back(vec2());
        edgesFaces0.push_back(nfaces);
        edgesFaces1.push_back(0);
        edgesNodes0.push_back(facesNodes1[nfaces]);
        edgesNodes1.push_back(facesNodes2[nfaces]);
        edgesLengths.push_back(0.);
        const uint nmin = std::min(facesNodes1[nfaces], facesNodes2[nfaces]);
        const uint nmax = std::max(facesNodes1[nfaces], facesNodes2[nfaces]);
        edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
    }
    if (!find_if_edge_in_mesh(facesNodes2[nfaces], facesNodes0[nfaces])) {
        edgesNormals.push_back(vec2());
        edgesCenters.push_back(vec2());
        edgesFaces0.push_back(nfaces);
        edgesFaces1.push_back(0);
        edgesNodes0.push_back(facesNodes2[nfaces]);
        edgesNodes1.push_back(facesNodes0[nfaces]);
        edgesLengths.push_back(0.);
        const uint nmin = std::min(facesNodes2[nfaces], facesNodes0[nfaces]);
        const uint nmax = std::max(facesNodes2[nfaces], facesNodes0[nfaces]);
        edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
    }
}



void mesh::read_su2_nodes(std::vector<std::string>& content) {
    auto nodes_line = find_tag("NPOIN", content);
    auto NPOIN = std::stoi(cut_str(content[nodes_line], '=')[1]);

    nodes.resize(NPOIN);
    for (int i=0; i<NPOIN; ++i) {
        auto np = str_to_floats( trim_all_whitespaces(content[nodes_line + 1 + i]) );
        nodes[i].x = np[0];
        nodes[i].y = np[1];
    }
}

void mesh::read_su2_elements(std::vector<std::string>& content) {
    auto start_line = find_tag("NELEM", content);
    auto NELEM = std::stoi(cut_str(content[start_line], '=')[1]);

    for (int i=0; i<NELEM; ++i) {
        auto np = str_to_ints( trim_all_whitespaces(content[start_line + 1 + i]) );

        // Add element
        facesAreas.push_back(0.);
        facesCenters.push_back(vec2());

        facesNodes0.push_back(np[1]);
        facesNodes1.push_back(np[2]);
        facesNodes2.push_back(np[3]);

        if (np[0] == 5) {
            // Triangle
            facesNodes3.push_back(0);
            faceIsTriangle.push_back(true);


            // Also add its edges
            add_last_tri_edges();
        } else if (np[0] == 9) {
            // Quad
            facesNodes3.push_back(np[4]);
            faceIsTriangle.push_back(false);

            // Also add its edges
            add_last_quad_edges();
        } else {
            throw std::invalid_argument( "invalid element number " + std::to_string(i) );
        }
    }
}

void mesh::read_su2_boundaries(std::vector<std::string>& content) {

    auto start_line = find_tag("NMARK", content);
    auto NMARK = std::stoi(cut_str(content[start_line], '=')[1]);

    int current_line = start_line + 1;

    for (int i=0; i<NMARK; ++i) {

        auto tag = trim(cut_str(content[current_line], '=')[1]);
        current_line += 1;
        auto NELM = std::stoi(cut_str(content[current_line], '=')[1]);
        current_line += 1;

        for (int j=0; j<NELM; ++j) {
            auto bound = str_to_ints( trim_all_whitespaces(content[current_line]) );

            boundsNormals.push_back(vec2());
            boundsCenters.push_back(vec2());
            boundsFaces.push_back(0);
            boundsLengths.push_back(0.);
            boundsNodes0.push_back(bound[1]);
            boundsNodes1.push_back(bound[2]);
            boundsTags.push_back(tag);
            const uint nmin = std::min(bound[1], bound[2]);
            const uint nmax = std::max(bound[1], bound[2]);
            boundsRef.insert({ std::make_tuple(nmin, nmax), boundsLengths.size()-1 });
            boundsFuncs.push_back(boundaryFuncMap().at("null"));
            
            current_line += 1;
        }
    }

}



void mesh::read_raw_su2(const std::string filename, const bool verb=false) {
	// Read su2 file into incomplete mesh data

	// Read file into array of strings
	std::vector<std::string> content;
	content = read(filename);

	if (content.size() != 0) {
		
        // Find ndime
        auto ndime_line = find_tag("NDIME", content);

        auto NDIME = std::stoi(cut_str(content[ndime_line], '=')[1]);

        if (verb) {
            std::cout << "NDIME = " << NDIME << std::endl;
        }

        if (NDIME == 2) {
            // OK, correct number of dimensions
            // read nodes
            read_su2_nodes(content);

            if (verb) {
                std::cout << "N NODES = " << nodes.size() << std::endl;
            }

            // read boundaries
            read_su2_boundaries(content);

            if (verb) {
                std::cout << "N BOUNDS = " << boundsNormals.size() << std::endl;
            }

            // read elements
            read_su2_elements(content);

            if (verb) {
                std::cout << "N ELEMS = " << faceIsTriangle.size() << std::endl;
            }


        } else {
            throw std::invalid_argument( "Mesh is not 2 dimensional" );
        }
	}
}


// Mesh constructor from file
mesh::mesh(const std::string filename, const bool verb=true) {
    // open file
    auto start = std::chrono::high_resolution_clock::now();

	if (verb) {
		std::cout << "(1/3) Reading mesh from file '" << filename << "' ... " << std::flush;
	}
	read_raw_su2(filename);
	if (verb) {
		std::cout << "File read. \n(2/3) Converting format to finite volume mesh ... " << std::flush;
	}
	// Convert node and face information
	convert_node_face_info();
	if (verb) {
		std::cout << "Conversion completed. \n(3/3) Computing mesh metrics ... " << std::flush;
	}
	// Compute normals, lengths and areas
	compute_mesh();
    auto stop = std::chrono::high_resolution_clock::now();

	if (verb) {
		std::cout << "Mesh import completed.\n" << std::flush;
	}

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    import_duration = ((double) duration.count())/1000.;

	if (verb) {

		std::cout << "\nMesh metrics\n";
		std::cout << " - Number of nodes = " << nodes.size() << "\n";
		std::cout << " - Number of edges = " << edgesNormals.size() << "\n";
		std::cout << " - Number of faces = " << facesAreas.size() << "\n";
        std::cout << " - Import time     = " << ((double) duration.count())/1000. << " s\n";
		std::cout << std::endl;
	}
}



#endif

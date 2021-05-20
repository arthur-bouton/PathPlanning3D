#include "path_search.h"

#include <queue>
#include <Eigen/Dense>
#include <chrono>


// Maximum number of extra entries that can be added to the node table:
#define EXTRA_VERTICES_MAX 3


PathSearch::PathSearch( const Eigen::MatrixXd& mesh_vertices,
                        const Eigen::MatrixXi& mesh_faces,
                        std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)> metric,
                        std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)> heuristic ) :
                        metric_( metric ),
                        heuristic_( heuristic )
{
	// Define the metric and heuristic functions used by default:
	auto default_metric = []( Eigen::RowVector3d vertex_1, Eigen::RowVector3d vertex_2 )
	{
		return ( vertex_1 - vertex_2 ).norm();
	};
	if ( !metric )
		metric_ = default_metric;
	if ( !heuristic )
		heuristic_ = default_metric;


	//---------------------------------------------//
	// Preprocess the mesh to build the node table //
	//---------------------------------------------//

	// Reserve enough memory for the node table to guarantee the validity of its pointers:
	node_table_.reserve( mesh_vertices.rows() + EXTRA_VERTICES_MAX );

	// Register every vertex of the mesh once:
	duplicate_table_ = std::vector<node_t*>( mesh_vertices.rows(), nullptr );
	for ( Eigen::Index vertex_id = 0 ; vertex_id < mesh_vertices.rows() ; ++vertex_id )
	{
		// Check if this vertex has already been registered:
		for ( node_t& node : node_table_ )
			if ( node.vertex == mesh_vertices.row( vertex_id ) )
			{
				// If it is a duplicate, register a pointer to its version stored in the node table:
				duplicate_table_[vertex_id] = &node;
				break;
			}

		if ( !duplicate_table_[vertex_id] )
		{
			// If it is a new vertex, add it to the tables:
			node_table_.push_back( node_t( mesh_vertices.row( vertex_id ) ) );
			duplicate_table_[vertex_id] = &node_table_.back();
		}
	}

	// Register the number of unique vertices in the mesh:
	n_mesh_vertices_ = node_table_.size();

	// Register the pointers to the adjacent nodes of each vertex:
	for ( Eigen::Index i = 0 ; i < mesh_faces.rows() ; ++i )
		for ( Eigen::Index j = 0 ; j < mesh_faces.cols() ; ++j )
		{
			// Original index of the considered vertex:
			unsigned int vertex_id = mesh_faces( i, j );
			// Pointer to the corresponding node in the node table:
			node_t* node_ptr = duplicate_table_[vertex_id];

			std::array<node_t*,2> face_nodes;

			for ( int k : { 0, 1 } )
			{
				// Original index of an adjacent vertex:
				unsigned int adjacent_vertex_id = mesh_faces( i, ( j + k + 1 )%3 );
				// Pointer to the corresponding node in the node table:
				node_t* adjacent_node_ptr = duplicate_table_[adjacent_vertex_id];
				// Add the adjacent-node pointer to the set of the current node:
				node_ptr->adjacent_nodes.insert( adjacent_node_ptr );
				// Store the pointer to the other node that compose this face:
				face_nodes[k] = adjacent_node_ptr;
			}
			// Add the face to the nodes that are part of it:
			node_ptr->faces.push_back( face_nodes );
		}
}


void PathSearch::reset_node_table()
{
	// Remove every link to the extra vertices that don't belong to the initial mesh:
	for ( std::size_t i = n_mesh_vertices_ ; i < node_table_.size() ; ++i )
		for ( node_t* adjacent_node_ptr : node_table_[i].adjacent_nodes )
			adjacent_node_ptr->adjacent_nodes.erase( &node_table_[i] );

	// Remove the extra vertices that don't belong to the initial mesh:
	node_table_.resize( n_mesh_vertices_, node_t( Eigen::RowVector3d() ) );

	// Reset each node:
	for ( node_t& node : node_table_ )
	{
		node.path.clear();
		node.status = node_t::unseen;
	}
}


PathSearch::node_t* PathSearch::findNode( const Eigen::RowVector3d& vertex )
{
	// Search for the vertex in the node table:
	for ( node_t& node : node_table_ )
		if ( node.vertex == vertex )
			return &node;

	// Check if we don't risk a reallocation of the node table:
	if ( node_table_.size() == node_table_.capacity() )
		throw std::runtime_error( "Too many vertices are different from the initial mesh vertices" );

	// If the vertex doesn't belong to the initial mesh, search which face it belongs to:
	node_t* face_nodes[3];
	findFaceNodes( vertex, face_nodes );
	
	// Add a temporary node to the node table corresponding to this extra vertex:
	node_t new_node( vertex );
	for ( node_t* face_node_ptr : face_nodes )
		new_node.adjacent_nodes.insert( face_node_ptr );
	node_table_.push_back( new_node );

	// Add the new node as an adjacent node for the nodes of the face it belongs to:
	for ( node_t* face_node_ptr : face_nodes )
		face_node_ptr->adjacent_nodes.insert( &node_table_.back() );

	return &node_table_.back();
}


void PathSearch::findFaceNodes( const Eigen::RowVector3d& vertex, PathSearch::node_t* face_nodes[3] )
{
	// Sort the mesh vertices by their distance from the given vertex:
	typedef struct dist_id { float distance; node_t* node_ptr; } dist_id_t;
	auto cmp = []( dist_id_t left, dist_id_t right ) { return left.distance > right.distance; };
	std::priority_queue<dist_id,std::vector<dist_id>,decltype(cmp)> distance_queue( cmp );

	for ( node_t& node : node_table_ )
	{
		dist_id_t new_entry;
		new_entry.node_ptr = &node;
		new_entry.distance = ( vertex - node.vertex ).norm();
		distance_queue.push( new_entry );
	}

	// Search among every node, starting by the closest ones to the vertex:
	while ( !distance_queue.empty() )
	{
		node_t* node_ptr = distance_queue.top().node_ptr;

		// Search among the faces the vertex of this node belongs to:
		for ( const std::array<node_t*,2>& face : node_ptr->faces )
		{
			// Gather the three vertices of this face:
			face_nodes[0] = node_ptr;
			face_nodes[1] = face[0];
			face_nodes[2] = face[1];

			// Compute the cross product with each edge of the current face:
			Eigen::Matrix3d normals;
			for ( int i = 0 ; i < 3 ; ++i )
				normals.row( i ) = ( face_nodes[i]->vertex - vertex ).cross( face_nodes[(i+1)%3]->vertex - face_nodes[i]->vertex );

			// Check if the cross-product results are pointing toward the same direction,
			// which means that the vertex is inside the current face
			// (if the cross products are also collinear, it means that the vertex lies exactly on the triangle surface):
			if ( normals.row( 0 ).dot( normals.row( 1 ) )*normals.row( 1 ).dot( normals.row( 2 ) ) >= 0
			  && normals.row( 1 ).dot( normals.row( 2 ) )*normals.row( 2 ).dot( normals.row( 3 ) ) >= 0 )
				return;
		}

		distance_queue.pop();
	}

	throw std::runtime_error( "The vertex don't belong to the mesh" );
}


Eigen::MatrixXd PathSearch::translatePath( const std::vector<node_t*>& path_nodes ) const
{
	// Fill the vertex matrix corresponding to the path:
	Eigen::MatrixXd path_vertices( path_nodes.size(), 3 );
	for ( std::size_t i = 0 ; i < path_nodes.size() ; ++i )
		path_vertices.row( i ) = path_nodes[i]->vertex;
	
	return path_vertices;
}


Eigen::MatrixXd PathSearch::A_star_search( const std::vector<Eigen::RowVector3d>& start_vertices,
                                           const Eigen::RowVector3d& goal_vertex )
{
	reset_node_table();

	// We use a multiset so that we don't have to scan the whole list to find the open node with the minimum f-score:
	auto cmp = [this]( node_t* left, node_t* right ) { return left->f_score < right->f_score; };
	// The nodes will be sorted according to their f-score:
	std::multiset<node_t*,decltype(cmp)> openNodes( cmp );

	for ( const Eigen::RowVector3d& start_vertex : start_vertices )
	{
		// Initialize the starting node:
		node_t* start_node_ptr = findNode( start_vertex );
		start_node_ptr->path.push_back( start_node_ptr );
		start_node_ptr->g_score = 0;
		start_node_ptr->f_score = heuristic_( start_vertex, goal_vertex );

		// Add the starting node to the queue of open nodes:
		openNodes.insert( start_node_ptr );
		start_node_ptr->status = node_t::open;
	}

	// Get the pointer to the goal node:
	node_t* goal_node_ptr = findNode( goal_vertex );


	// Initialize the clock to watch the performance:
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	// Total number of nodes opened:
	unsigned int n_opened = 0;


	while ( !openNodes.empty() )
	{
		// Take the node with the minimum f-score among the open nodes:
		node_t* current_node_ptr = *openNodes.begin();
		++n_opened;

		// If we are expanding the target node, it means that we have now the optimal path:
		if ( current_node_ptr == goal_node_ptr )
		{
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			printf( "Target reached in %6ld µs | nodes opened: %5u | path length: %3zu vertices | path cost: %f\n",
					std::chrono::duration_cast<std::chrono::microseconds>( end - begin ).count(),
					++n_opened, current_node_ptr->path.size(), current_node_ptr->g_score );
			return translatePath( current_node_ptr->path );
		}

		// Remove the current node from the queue of open nodes:
		openNodes.erase( openNodes.begin() );
		current_node_ptr->status = node_t::closed;

		// Expand the adjacent nodes:
		for ( node_t* adjacent_node_ptr : current_node_ptr->adjacent_nodes )
		{
			// Compute the backward cost of this adjacent node:
			float new_g_score = current_node_ptr->g_score + metric_( current_node_ptr->vertex, adjacent_node_ptr->vertex );

			// Check if this adjacent node has been seen before:
			if ( adjacent_node_ptr->status != node_t::unseen )
			{
				// Check if the new path to this node is cheaper than previously:
				// Always true if the node has already been closed and the heuristic is consistent.
				// Otherwise, reopening it ensures the optimality of the final path.
				if ( new_g_score > adjacent_node_ptr->g_score )
					continue;
				// Check if the node has been marked as open to potentially avoid the cost of searching for it in openNodes:
				else if ( adjacent_node_ptr->status == node_t::open )
					// Keep only the current best path to this adjacent node:
					openNodes.erase( adjacent_node_ptr );
			}

			adjacent_node_ptr->path = current_node_ptr->path;
			adjacent_node_ptr->path.push_back( adjacent_node_ptr );
			adjacent_node_ptr->g_score = new_g_score;
			adjacent_node_ptr->f_score = new_g_score + heuristic_( adjacent_node_ptr->vertex, goal_vertex );

			// Add this node to the queue of nodes to expand:
			openNodes.insert( adjacent_node_ptr );
			adjacent_node_ptr->status = node_t::open;
		}
	}

	// No path to the target has been found:
	return goal_vertex;
}


Eigen::VectorXd PathSearch::dijkstra_scan( const std::vector<Eigen::RowVector3d>& sources )
{
	reset_node_table();

	// We use a multiset so that we don't have to scan the whole list to find the open node with the minimum g-score:
	auto cmp = [this]( node_t* left, node_t* right ) { return left->g_score < right->g_score; };
	// The nodes will be sorted according to their g-score:
	std::multiset<node_t*,decltype(cmp)> openNodes( cmp );


	for ( const Eigen::RowVector3d& source_vertex : sources )
	{
		// Initialize the source node:
		node_t* source_node_ptr = findNode( source_vertex );
		source_node_ptr->g_score = 0;

		// Add the source node to the queue of open nodes:
		openNodes.insert( source_node_ptr );
		source_node_ptr->status = node_t::open;
	}


	while ( !openNodes.empty() )
	{
		// Take the node with the minimum g-score among the open nodes:
		node_t* current_node_ptr = *openNodes.begin();

		// Remove the current node from the queue of open nodes:
		openNodes.erase( openNodes.begin() );
		current_node_ptr->status = node_t::closed;

		// Expand the adjacent nodes:
		for ( node_t* adjacent_node_ptr : current_node_ptr->adjacent_nodes )
		{
			// Compute the backward cost of this adjacent node:
			float new_g_score = current_node_ptr->g_score + metric_( current_node_ptr->vertex, adjacent_node_ptr->vertex );

			// Check if this adjacent node has been seen before:
			if ( adjacent_node_ptr->status != node_t::unseen )
			{
				// Check if the new path to this node is cheaper than previously:
				if ( new_g_score > adjacent_node_ptr->g_score )
					continue;
				// Check if the node has been marked as open to potentially avoid the cost of searching for it in openNodes:
				else if ( adjacent_node_ptr->status == node_t::open )
					// Keep only the current shortest path to this adjacent node:
					openNodes.erase( adjacent_node_ptr );
			}

			adjacent_node_ptr->g_score = new_g_score;

			// Add this node to the queue of nodes to expand:
			openNodes.insert( adjacent_node_ptr );
			adjacent_node_ptr->status = node_t::open;
		}
	}

	// Fill a distance vector using the original vertex indices:
	Eigen::VectorXd distances( duplicate_table_.size() );
	for ( std::size_t i = 0 ; i < duplicate_table_.size() ; ++i )
		distances[i] = duplicate_table_[i]->g_score;

	return distances;
}

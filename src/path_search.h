#ifndef PATH_SEARCH_HH
#define PATH_SEARCH_HH 

#include <Eigen/Core>
#include <vector>
#include <set>
#include <array>
#include <memory>


class PathSearch
{
	public:

	typedef std::shared_ptr<PathSearch> ptr_t;

	PathSearch( const Eigen::MatrixXd& mesh_vertices,
	            const Eigen::MatrixXi& mesh_faces,
	            std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)> metric = nullptr,
	            std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)> heuristic = nullptr );

	// Return the shortest path as a matrix of vertices from the closest starting vertex of start_vertices to goal_vertex:
	Eigen::MatrixXd A_star_search( const std::vector<Eigen::RowVector3d>& start_vertices,
	                               const Eigen::RowVector3d& goal_vertex );

	// Return the distance of each vertex from the source(s):
	Eigen::VectorXd dijkstra_scan( const std::vector<Eigen::RowVector3d>& sources );

	inline void setMetric(    const std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)>& metric )    { metric_ = metric; };
	inline void setHeuristic( const std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)>& heuristic ) { heuristic_ = heuristic; };

	protected:

	typedef struct node
	{
		
		Eigen::RowVector3d vertex;                     // Corresponding vertex.
		std::set<struct node*> adjacent_nodes;         // Pointers to the adjacent nodes. The set ensures their uniqueness.
		std::vector<struct node*> path;                // Current best path to this node.
		float g_score, f_score;                        // Backward cost and priority of the node.
		enum { unseen, open, closed } status = unseen; // Current status of the node.
		std::vector<std::array<struct node*,2>> faces; // Pointers the nodes of each face this node is part of
		                                               // (only used when a vertex is not part of the mesh).
		node( Eigen::RowVector3d vert ) : vertex( vert ) {}
	} node_t;

	// Table that stores all the nodes, each corresponding to a unique vertex of the mesh:
	std::vector<node_t> node_table_;

	void reset_node_table();

	node_t* findNode( const Eigen::RowVector3d& vertex );
	void findFaceNodes( const Eigen::RowVector3d& vertex, node_t* face_nodes[3] );
	Eigen::MatrixXd translatePath( const std::vector<node_t*>& path_nodes ) const;

	std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)> metric_;
	std::function<float(Eigen::MatrixXd,Eigen::MatrixXd)> heuristic_;

	// Table that translates the original vertex indices to the corresponding node pointers:
	std::vector<node_t*> duplicate_table_;
	// Number of unique vertices in the initial the mesh:
	size_t n_mesh_vertices_ = 0;
};


#endif

#define MAX_ELEMENT_VERTS 4

enum FV_MESH_GMESH_ELEMENT
{
	FV_MESH_MIN_GMESH_ELEMENT,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_LINE = 1,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_TRIANGLE = 2,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_QUADRANGLE = 3,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_TETRAHEDRAL = 4,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_HEXAHEDRAL = 5,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_PRISM = 6,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_PYRAMID = 7,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_LINE = 8,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_TRIANGLE = 9,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_QUADRANGLE = 10,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_TETRAHEDRAL = 11,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_HEXAHEDRAL = 12,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_PRISM = 13,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_PYRAMID = 14,
	FV_MESH_GMESH_ELEMENT_POINT = 15,
	FV_MESH_MAX_GMESH_ELEMENT
};


//using namespace std;


typedef struct
{
    double r[3];
    int mid;
}
vertex_t;

typedef struct
{
    int vertex_mids[MAX_ELEMENT_VERTS];
    int nvert;
    int mid;
    int type;
}
element_t;

class Mesh
{
public:
    Mesh(std::string filename);
    ~Mesh();

    void writeMesh(std::string filename, std::string value_label, const vector<double>& values) const;
    int getElement3DCount() const {return elements_3d_.size(); }

    //void partitionMesh();                           //write this function
    //void outputStatistics() const;                  //write this function

private:

    void readMesh(std::string filename);
    void getElementVertices();
    const vertex_t* getVertexFromMID(int vertex_mid) const;
    void createElementCentroidsList();


    int rank_;
    int nproc_;
    vector<vertex_t> vertices_;
    vector<element_t> elements_3d_;
    std::vector<double3_t> element_3d_centroids_;
};

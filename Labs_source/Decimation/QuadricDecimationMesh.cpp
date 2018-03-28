/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode
    QuadricDecimationMesh::QuadricIsoSurfaces =
        NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize() {
  // Allocate memory for the quadric array
  size_t numVerts = mVerts.size();
  mQuadrics.reserve(numVerts);
  std::streamsize width = std::cerr.precision(); // store stream precision
  for (size_t i = 0; i < numVerts; i++) {

    // Compute quadric for vertex i here
    mQuadrics.push_back(createQuadricForVert(i));

    // Calculate initial error, should be numerically close to 0

    Vector3<float> v0 = mVerts[i].pos;
    Vector4<float> v(v0[0], v0[1], v0[2], 1);
    Matrix4x4<float> m = mQuadrics.back();

    float error = v * (m * v);
    // std::cerr << std::scientific << std::setprecision(2) << error << " ";
  }
  std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

  // Run the initialize for the parent class to initialize the edge collapses
  DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute,
 * DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse *collapse) {
  // Compute collapse->position and collapse->cost here
  // based on the quadrics at the edge endpoints
	size_t idx1 = e(collapse->halfEdge).vert;
	size_t idx2 = e(e(collapse->halfEdge).pair).vert;
	Matrix4x4<float> Q1 = createQuadricForVert(idx1);
	Matrix4x4<float> Q2 = createQuadricForVert(idx2);
	Matrix4x4<float> Q = Q1 + Q2;
	Matrix4x4<float> Qpos = Q;
	Qpos(3, 0) = Qpos(3,1) = Qpos(3,2) = 0.0f;
	Qpos(3, 3) = 1.0f;

	if (!Qpos.IsSingular())
	{
		//Find the new minimum collapse position
		Vector4<float> homoV = { 0.0f, 0.0f , 0.0f , 1.0f };
		Vector4<float> V;
		V = Qpos.Inverse()*homoV;
		collapse->position = Vector3<float>(V[0], V[1], V[2]);

		//Calculate the cost
		float C = V * (Q * V);
		collapse->cost = C;

		
	}
	else
	{
		Vector3<float> v1 = v(idx1).pos;
		Vector3<float> v2 = v(idx2).pos;
		Vector3<float> midPoint = (v1 + v2) / 2.0f;

		Vector4<float> v1x4 = { v1[0], v1[1], v1[2], 1.0f };
		Vector4<float> v2x4 = { v2[0], v2[1], v2[2], 1.0f };
		Vector4<float> vMidx4 = { midPoint[0], midPoint[1], midPoint[2], 1.0f };

		float c1 = v1x4 * (Q*v1x4);
		float c2 = v2x4 * (Q*v2x4);
		float cMid = vMidx4 * (Q*vMidx4);

		collapse->cost = std::min(std::min(c1, c2), cMid);

		if (collapse->cost == c1)
			collapse->position = v1;

		else if (collapse->cost == c2)
			collapse->position = v2;

		else
			collapse->position = midPoint;
	}
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(size_t ind) {
  DecimationMesh::updateVertexProperties(ind);
  mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
Matrix4x4<float>
QuadricDecimationMesh::createQuadricForVert(size_t indx) const {
  float q[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  Matrix4x4<float> Q(q);
 
  // The quadric for a vertex is the sum of all the quadrics for the adjacent
  // faces Tip: Matrix4x4 has an operator +=
  std::vector<size_t> foundFaces = FindNeighborFaces(indx);
  for (size_t i = 0; i < foundFaces.size(); i++)
  {
	  Q += createQuadricForFace(foundFaces.at(i));
  }
  return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
Matrix4x4<float>
QuadricDecimationMesh::createQuadricForFace(size_t indx) const {

  // Calculate the quadric (outer product of plane parameters) for a face
  // here using the formula from Garland and Heckbert

	Vector3<float> p = v(e(f(indx).edge).vert).pos;
	Vector3<float> n = f(indx).normal;

	float d = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
	d = -1.0*d;

	float q[4][4] = { 
		{ n[0] * n[0], n[0] * n[1], n[0] * n[2], n[0] * d },
		{ n[1] * n[0], n[1] * n[1], n[1] * n[2], n[1] * d },
		{ n[2] * n[0], n[2] * n[1], n[2] * n[2], n[2] * d },
		{ n[0] * d, n[1] * d , n[2] * d, d * d }
					
	};
  return Matrix4x4<float>(q);
}

void QuadricDecimationMesh::Render() {
  DecimationMesh::Render();

  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  if (mVisualizationMode == QuadricIsoSurfaces) {
    // Apply transform
    glPushMatrix(); // Push modelview matrix onto stack

    // Implement the quadric visualization here
    std::cout << "Quadric visualization not implemented" << std::endl;

    // Restore modelview matrix
    glPopMatrix();
  }
}

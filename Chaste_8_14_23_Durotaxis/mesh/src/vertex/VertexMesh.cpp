/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "VertexMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements)
        : mpDelaunayMesh(nullptr)
{

    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements
    for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index = 0; elem_index < vertexElements.size(); elem_index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // In 3D, populate mFaces
    if (SPACE_DIM == 3)
    {
        // Use a std::set to keep track of which faces have been added to mFaces
        std::set<unsigned> faces_counted;

        // Loop over mElements
        for (unsigned elem_index = 0; elem_index < mElements.size(); elem_index++)
        {
            // Loop over faces of this element
            for (unsigned face_index = 0; face_index < mElements[elem_index]->GetNumFaces(); face_index++)
            {
                VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = mElements[elem_index]->GetFace(face_index);
                unsigned global_index = p_face->GetIndex();

                // If this face is not already contained in mFaces, add it, and update faces_counted
                if (faces_counted.find(global_index) == faces_counted.end())
                {
                    mFaces.push_back(p_face);
                    faces_counted.insert(global_index);
                }
            }
        }
    }

    // Register elements with nodes
    for (unsigned index = 0; index < mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[index];

        unsigned element_index = p_element->GetIndex();
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned node_index = 0; node_index < num_nodes_in_element; node_index++)
        {
            p_element->GetNode(node_index)->AddElement(element_index);
        }
    }

    this->mMeshChangesDuringSimulation = false;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                               std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces,
                                               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements)
        : mpDelaunayMesh(nullptr)
{
    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    // Populate mNodes mFaces and mElements
    for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }

    for (unsigned face_index = 0; face_index < faces.size(); face_index++)
    {
        VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_temp_face = faces[face_index];
        mFaces.push_back(p_temp_face);
    }

    for (unsigned elem_index = 0; elem_index < vertexElements.size(); elem_index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // Register elements with nodes
    for (unsigned index = 0; index < mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index = 0; node_index < p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = false;
}

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
VertexMesh<2, 2>::VertexMesh(TetrahedralMesh<2, 2>& rMesh, bool isPeriodic)
        : mpDelaunayMesh(&rMesh)
{
    //Note  !isPeriodic is not used except through polymorphic calls in rMesh

    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    unsigned num_elements = mpDelaunayMesh->GetNumAllNodes();
    unsigned num_nodes = mpDelaunayMesh->GetNumAllElements();

    // Allocate memory for mNodes and mElements
    this->mNodes.reserve(num_nodes);

    // Create as many elements as there are nodes in the mesh
    mElements.reserve(num_elements);
    for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
    {
        VertexElement<2, 2>* p_element = new VertexElement<2, 2>(elem_index);
        mElements.push_back(p_element);
    }

    // Populate mNodes
    GenerateVerticesFromElementCircumcentres(rMesh);

    // Loop over elements of the Delaunay mesh (which are nodes/vertices of this mesh)
    for (unsigned i = 0; i < num_nodes; i++)
    {
        // Loop over nodes owned by this triangular element in the Delaunay mesh
        // Add this node/vertex to each of the 3 vertex elements
        for (unsigned local_index = 0; local_index < 3; local_index++)
        {
            unsigned elem_index = mpDelaunayMesh->GetElement(i)->GetNodeGlobalIndex(local_index);
            unsigned num_nodes_in_elem = mElements[elem_index]->GetNumNodes();
            unsigned end_index = num_nodes_in_elem > 0 ? num_nodes_in_elem - 1 : 0;

            mElements[elem_index]->AddNode(this->mNodes[i], end_index);
        }
    }

    // Reorder mNodes anticlockwise
    for (unsigned elem_index = 0; elem_index < mElements.size(); elem_index++)
    {
        /**
         * Create a std::vector of pairs, where each pair comprises the angle
         * between the centre of the Voronoi element and each node with that
         * node's global index in the Voronoi mesh.
         */
        std::vector<std::pair<double, unsigned> > index_angle_list;
        for (unsigned local_index = 0; local_index < mElements[elem_index]->GetNumNodes(); local_index++)
        {
            c_vector<double, 2> vectorA = mpDelaunayMesh->GetNode(elem_index)->rGetLocation();
            c_vector<double, 2> vectorB = mElements[elem_index]->GetNodeLocation(local_index);
            c_vector<double, 2> centre_to_vertex = mpDelaunayMesh->GetVectorFromAtoB(vectorA, vectorB);

            double angle = atan2(centre_to_vertex(1), centre_to_vertex(0));
            unsigned global_index = mElements[elem_index]->GetNodeGlobalIndex(local_index);

            std::pair<double, unsigned> pair(angle, global_index);
            index_angle_list.push_back(pair);
        }

        // Sort the list in order of increasing angle
        sort(index_angle_list.begin(), index_angle_list.end());

        // Create a new Voronoi element and pass in the appropriate Nodes, ordered anticlockwise
        VertexElement<2, 2>* p_new_element = new VertexElement<2, 2>(elem_index);
        for (unsigned count = 0; count < index_angle_list.size(); count++)
        {
            unsigned local_index = count > 1 ? count - 1 : 0;
            p_new_element->AddNode(mNodes[index_angle_list[count].second], local_index);
        }

        // Replace the relevant member of mElements with this Voronoi element
        delete mElements[elem_index];
        mElements[elem_index] = p_new_element;
    }

    this->mMeshChangesDuringSimulation = false;
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
VertexMesh<3, 3>::VertexMesh(TetrahedralMesh<3, 3>& rMesh)
        : mpDelaunayMesh(&rMesh)
{
    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    unsigned num_nodes = mpDelaunayMesh->GetNumAllElements();

    // Allocate memory for mNodes
    this->mNodes.reserve(num_nodes);

    // Populate mNodes
    GenerateVerticesFromElementCircumcentres(rMesh);

    std::map<unsigned, VertexElement<3, 3>*> index_element_map;
    unsigned face_index = 0;
    unsigned element_index = 0;

    // Loop over each edge of the Delaunay mesh and populate mFaces and mElements
    for (TetrahedralMesh<3, 3>::EdgeIterator edge_iterator = mpDelaunayMesh->EdgesBegin();
         edge_iterator != mpDelaunayMesh->EdgesEnd();
         ++edge_iterator)
    {
        Node<3>* p_node_a = edge_iterator.GetNodeA();
        Node<3>* p_node_b = edge_iterator.GetNodeB();

        if (!(p_node_a->IsBoundaryNode() && p_node_b->IsBoundaryNode()))
        {
            std::set<unsigned>& node_a_element_indices = p_node_a->rGetContainingElementIndices();
            std::set<unsigned>& node_b_element_indices = p_node_b->rGetContainingElementIndices();
            std::set<unsigned> edge_element_indices;

            std::set_intersection(node_a_element_indices.begin(),
                                  node_a_element_indices.end(),
                                  node_b_element_indices.begin(),
                                  node_b_element_indices.end(),
                                  std::inserter(edge_element_indices, edge_element_indices.begin()));

            c_vector<double, 3> edge_vector;
            edge_vector = p_node_b->rGetLocation() - p_node_a->rGetLocation();

            c_vector<double, 3> mid_edge;
            mid_edge = edge_vector * 0.5 + p_node_a->rGetLocation();

            unsigned element0_index = *(edge_element_indices.begin());

            c_vector<double, 3> basis_vector1;
            basis_vector1 = mNodes[element0_index]->rGetLocation() - mid_edge;

            c_vector<double, 3> basis_vector2;
            basis_vector2 = VectorProduct(edge_vector, basis_vector1);

            /**
             * Create a std::vector of pairs, where each pair comprises the angle
             * between the centre of the Voronoi element and each node with that
             * node's global index in the Voronoi mesh.
             */
            std::vector<std::pair<double, unsigned> > index_angle_list;

            // Loop over each element containing this edge (i.e. those containing both nodes of the edge)
            for (std::set<unsigned>::iterator index_iter = edge_element_indices.begin();
                 index_iter != edge_element_indices.end();
                 ++index_iter)
            {
                // Calculate angle
                c_vector<double, 3> vertex_vector = mNodes[*index_iter]->rGetLocation() - mid_edge;

                double local_vertex_dot_basis_vector1 = inner_prod(vertex_vector, basis_vector1);
                double local_vertex_dot_basis_vector2 = inner_prod(vertex_vector, basis_vector2);

                double angle = atan2(local_vertex_dot_basis_vector2, local_vertex_dot_basis_vector1);

                std::pair<double, unsigned> pair(angle, *index_iter);
                index_angle_list.push_back(pair);
            }

            // Sort the list in order of increasing angle
            sort(index_angle_list.begin(), index_angle_list.end());

            // Create face
            VertexElement<2, 3>* p_face = new VertexElement<2, 3>(face_index);
            face_index++;
            for (unsigned count = 0; count < index_angle_list.size(); count++)
            {
                unsigned local_index = count > 1 ? count - 1 : 0;
                p_face->AddNode(mNodes[index_angle_list[count].second], local_index);
            }

            // Add face to list of faces
            mFaces.push_back(p_face);

            // Add face to appropriate elements
            if (!p_node_a->IsBoundaryNode())
            {
                unsigned node_a_index = p_node_a->GetIndex();
                if (index_element_map[node_a_index])
                {
                    // If there is already an element with this index, add the face to it...
                    index_element_map[node_a_index]->AddFace(p_face);
                }
                else
                {
                    // ...otherwise create an element, add the face to it, and add to the map
                    mVoronoiElementIndexMap[node_a_index] = element_index;
                    VertexElement<3, 3>* p_element = new VertexElement<3, 3>(element_index);
                    element_index++;
                    p_element->AddFace(p_face);
                    index_element_map[node_a_index] = p_element;
                }
            }
            if (!p_node_b->IsBoundaryNode())
            {
                unsigned node_b_index = p_node_b->GetIndex();
                if (index_element_map[node_b_index])
                {
                    // If there is already an element with this index, add the face to it...
                    index_element_map[node_b_index]->AddFace(p_face);
                }
                else
                {
                    // ...otherwise create an element, add the face to it, and add to the map
                    mVoronoiElementIndexMap[node_b_index] = element_index;
                    VertexElement<3, 3>* p_element = new VertexElement<3, 3>(element_index);
                    element_index++;
                    p_element->AddFace(p_face);
                    index_element_map[node_b_index] = p_element;
                }
            }
        }
    }

    // Populate mElements
    unsigned elem_count = 0;
    for (std::map<unsigned, VertexElement<3, 3>*>::iterator element_iter = index_element_map.begin();
         element_iter != index_element_map.end();
         ++element_iter)
    {
        mElements.push_back(element_iter->second);
        mVoronoiElementIndexMap[element_iter->first] = elem_count;
        elem_count++;
    }

    this->mMeshChangesDuringSimulation = false;
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::GenerateVerticesFromElementCircumcentres(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    c_matrix<double, ELEMENT_DIM, SPACE_DIM> inverse_jacobian;
    double jacobian_det;

    // Loop over elements of the Delaunay mesh and populate mNodes
    for (unsigned i = 0; i < rMesh.GetNumElements(); i++)
    {
        // Calculate the circumcentre of this element in the Delaunay mesh
        rMesh.GetInverseJacobianForElement(i, jacobian, jacobian_det, inverse_jacobian);
        c_vector<double, SPACE_DIM + 1> circumsphere = rMesh.GetElement(i)->CalculateCircumsphere(jacobian, inverse_jacobian);

        c_vector<double, SPACE_DIM> circumcentre;
        for (unsigned j = 0; j < SPACE_DIM; j++)
        {
            circumcentre(j) = circumsphere(j);
        }

        // Create a node in the Voronoi mesh at the location of this circumcentre
        this->mNodes.push_back(new Node<SPACE_DIM>(i, circumcentre));
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetEdgeLength(unsigned elementIndex1, unsigned elementIndex2)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    std::set<unsigned> node_indices_1;
    for (unsigned i = 0; i < mElements[elementIndex1]->GetNumNodes(); i++)
    {
        node_indices_1.insert(mElements[elementIndex1]->GetNodeGlobalIndex(i));
    }
    std::set<unsigned> node_indices_2;
    for (unsigned i = 0; i < mElements[elementIndex2]->GetNumNodes(); i++)
    {
        node_indices_2.insert(mElements[elementIndex2]->GetNodeGlobalIndex(i));
    }

    std::set<unsigned> shared_nodes;
    std::set_intersection(node_indices_1.begin(), node_indices_1.end(),
                          node_indices_2.begin(), node_indices_2.end(),
                          std::inserter(shared_nodes, shared_nodes.begin()));

    if (shared_nodes.size() == 1)
    {
        // It's possible that these two elements are actually infinite but are on the edge of the domain
        EXCEPTION("Elements " << elementIndex1 << " and " << elementIndex2 << " share only one node.");
    }
    assert(shared_nodes.size() == 2);

    unsigned index1 = *(shared_nodes.begin());
    unsigned index2 = *(++(shared_nodes.begin()));

    double edge_length = this->GetDistanceBetweenNodes(index1, index2);
    return edge_length;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElongationShapeFactorOfElement(unsigned index)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    double discriminant = sqrt((moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2));

    // Note that as the matrix of second moments of area is symmetric, both its eigenvalues are real
    double largest_eigenvalue = (moments(0) + moments(1) + discriminant) * 0.5;
    double smallest_eigenvalue = (moments(0) + moments(1) - discriminant) * 0.5;

    double elongation_shape_factor = sqrt(largest_eigenvalue / smallest_eigenvalue);
    return elongation_shape_factor;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh()
{
    mpDelaunayMesh = nullptr;
    this->mMeshChangesDuringSimulation = false;
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::~VertexMesh()
{
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
	//std::cout<<this->mNodes.size()<<"\t"<<index<<std::endl;
	//std::cout<<this<<"\t"<<index<<std::endl;
    assert(index < this->mNodes.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    ///\todo sort out boundary elements in a vertex mesh (#1263)
    //    assert(index < this->mBoundaryElements.size() );
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(unsigned elementIndex)
{
    unsigned node_index = UNSIGNED_UNSET;

    if (mVoronoiElementIndexMap.empty())
    {
        node_index = elementIndex;
    }
    else
    {
        for (std::map<unsigned, unsigned>::iterator iter = mVoronoiElementIndexMap.begin();
             iter != mVoronoiElementIndexMap.end();
             ++iter)
        {
            if (iter->second == elementIndex)
            {
                node_index = iter->first;
                break;
            }
        }
    }
    assert(node_index != UNSIGNED_UNSET);
    return node_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(unsigned nodeIndex)
{
    unsigned element_index = UNSIGNED_UNSET;

    if (mVoronoiElementIndexMap.empty())
    {
        element_index = nodeIndex;
    }
    else
    {
        std::map<unsigned, unsigned>::iterator iter = mVoronoiElementIndexMap.find(nodeIndex);

        if (iter == mVoronoiElementIndexMap.end())
        {
            EXCEPTION("This index does not correspond to a VertexElement");
        }
        else
        {
            element_index = iter->second;
        }
    }
    assert(element_index != UNSIGNED_UNSET);
    return element_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetRosetteRankOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    // Loop over nodes in the current element and find which is contained in the most elements
    unsigned rosette_rank = 0;
    for (unsigned node_idx = 0; node_idx < p_element->GetNumNodes(); node_idx++)
    {
        unsigned num_elems_this_node = p_element->GetNode(node_idx)->rGetContainingElementIndices().size();

        if (num_elems_this_node > rosette_rank)
        {
            rosette_rank = num_elems_this_node;
        }
    }

    // Return the rosette rank
    return rosette_rank;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // Delete elements
    for (unsigned i = 0; i < mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete faces
    for (unsigned i = 0; i < mFaces.size(); i++)
    {
        delete mFaces[i];
    }
    mFaces.clear();

    // Delete nodes
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    //std::cout<<"index\t"<<index<<"\t"<<"mElements.size()\t"<<mElements.size()<<std::endl;
    assert(index < mElements.size());
    return mElements[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
    return mFaces[index];
}


//added 9-21-21
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLengthOfElement_ab(unsigned index)
{
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    //unsigned num_nodes = p_element->GetNumNodes();

    //double length = zero_vector<double>(SPACE_DIM);
	c_vector<double, SPACE_DIM> length = zero_vector<double>(SPACE_DIM);

	length =  (p_element->GetNodeLocation(0) - p_element->GetNodeLocation(1));

	return abs(length(0));

}



//added 8-23-21
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfElement_ab(unsigned index)
{
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    //unsigned num_nodes = p_element->GetNumNodes();

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);


	centroid = 0.5 * (p_element->GetNodeLocation(0) + p_element->GetNodeLocation(1));

	return centroid;

}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfElement(unsigned index)
{
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    unsigned num_nodes = p_element->GetNumNodes();

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

    switch (SPACE_DIM)
    {
        case 1:
        {
            centroid = 0.5 * (p_element->GetNodeLocation(0) + p_element->GetNodeLocation(1));
        }
        break;
        case 2:
        {
            double centroid_x = 0;
            double centroid_y = 0;

            // Note that we cannot use GetVolumeOfElement() below as it returns the absolute, rather than signed, area
            double element_signed_area = 0.0;

            // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
            c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
            c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

            // Loop over vertices
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
                c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

                double this_x = pos_1[0];
                double this_y = pos_1[1];
                double next_x = pos_2[0];
                double next_y = pos_2[1];

                double signed_area_term = this_x * next_y - this_y * next_x;

                centroid_x += (this_x + next_x) * signed_area_term;
                centroid_y += (this_y + next_y) * signed_area_term;
                element_signed_area += 0.5 * signed_area_term;

                pos_1 = pos_2;
            }

            assert(element_signed_area != 0.0);

            // Finally, map back and employ GetVectorFromAtoB() to allow for periodicity
            centroid = first_node_location;
            centroid(0) += centroid_x / (6.0 * element_signed_area);
            centroid(1) += centroid_y / (6.0 * element_signed_area);
        }
        break;
        case 3:
        {
            ///\todo compute centroid rather than centre of mass (see #1422)
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                centroid += p_element->GetNodeLocation(local_index);
            }
            centroid /= ((double)num_nodes);
        }
        break;
        default:
            NEVER_REACHED;
    }
    return centroid;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = this->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Find the local index of this node in this element
        unsigned local_index = GetElement(*elem_iter)->GetNodeLocalIndex(nodeIndex);

        // Find the global indices of the preceding and successive nodes in this element
        unsigned num_nodes = GetElement(*elem_iter)->GetNumNodes();
        unsigned previous_local_index = (local_index + num_nodes - 1) % num_nodes;
        unsigned next_local_index = (local_index + 1) % num_nodes;

        // Add the global indices of these two nodes to the set of neighbouring node indices
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(previous_local_index));
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(next_local_index));
    }

    return neighbouring_node_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex)
{
    // Get a pointer to this element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(elemIndex);

    // Get the indices of the nodes contained in this element
    std::set<unsigned> node_indices_in_this_element;
    for (unsigned local_index = 0; local_index < p_element->GetNumNodes(); local_index++)
    {
        unsigned global_index = p_element->GetNodeGlobalIndex(local_index);
        node_indices_in_this_element.insert(global_index);
    }

    // Check that the node is in fact contained in the element
    if (node_indices_in_this_element.find(nodeIndex) == node_indices_in_this_element.end())
    {
        EXCEPTION("The given node is not contained in the given element.");
    }

    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices_not_in_this_element;

    // Get the indices of this node's neighbours
    std::set<unsigned> node_neighbours = GetNeighbouringNodeIndices(nodeIndex);

    // Check if each neighbour is also in this element; if not, add it to the set
    for (std::set<unsigned>::iterator iter = node_neighbours.begin();
         iter != node_neighbours.end();
         ++iter)
    {
        if (node_indices_in_this_element.find(*iter) == node_indices_in_this_element.end())
        {
            neighbouring_node_indices_not_in_this_element.insert(*iter);
        }
    }

    return neighbouring_node_indices_not_in_this_element;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringElementIndices(unsigned elementIndex)
{
    // Get a pointer to this element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(elementIndex);

    // Create a set of neighbouring element indices
    std::set<unsigned> neighbouring_element_indices;

    // Loop over nodes owned by this element
    for (unsigned local_index = 0; local_index < p_element->GetNumNodes(); local_index++)
    {
        // Get a pointer to this node
        Node<SPACE_DIM>* p_node = p_element->GetNode(local_index);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_node->rGetContainingElementIndices();

        // Form the union of this set with the current set of neighbouring element indices
        std::set<unsigned> all_elements;
        std::set_union(neighbouring_element_indices.begin(), neighbouring_element_indices.end(),
                       containing_elem_indices.begin(), containing_elem_indices.end(),
                       std::inserter(all_elements, all_elements.begin()));

        // Update the set of neighbouring element indices
        neighbouring_element_indices = all_elements;
    }

    // Lastly remove this element's index from the set of neighbouring element indices
    neighbouring_element_indices.erase(elementIndex);

    return neighbouring_element_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetMeshForVtk()
{
    return this;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void VertexMesh<1, 1>::ConstructFromMeshReader(AbstractMeshReader<1, 1>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("VertexMesh<1,1>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void VertexMesh<1, 2>::ConstructFromMeshReader(AbstractMeshReader<1, 2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("VertexMesh<1,2>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void VertexMesh<1, 3>::ConstructFromMeshReader(AbstractMeshReader<1, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("VertexMesh<1,3>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void VertexMesh<2, 3>::ConstructFromMeshReader(AbstractMeshReader<2, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("VertexMesh<2,3>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void VertexMesh<2, 2>::ConstructFromMeshReader(AbstractMeshReader<2, 2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(rMeshReader.HasNodePermutation() == false);
    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i = 0; i < num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool)node_data[2];
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for elements
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j = 0; j < num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        VertexElement<2, 2>* p_element = new VertexElement<2, 2>(elem_index, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = (unsigned)element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void VertexMesh<3, 3>::ConstructFromMeshReader(AbstractMeshReader<3, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(rMeshReader.HasNodePermutation() == false);

    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i = 0; i < num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool)node_data[3];
        node_data.pop_back();
        this->mNodes.push_back(new Node<3>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rMeshReader.GetNumElements());

    // Use a std::set to keep track of which faces have been added to mFaces
    std::set<unsigned> faces_counted;

    // Add elements
    for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
    {
        ///\todo Horrible hack! (#1076/#1377)
        typedef VertexMeshReader<3, 3> VERTEX_MESH_READER;
        assert(dynamic_cast<VERTEX_MESH_READER*>(&rMeshReader) != nullptr);

        // Get the data for this element
        VertexElementData element_data = static_cast<VERTEX_MESH_READER*>(&rMeshReader)->GetNextElementDataWithFaces();

        // Get the nodes owned by this element
        std::vector<Node<3>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j = 0; j < num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Get the faces owned by this element
        std::vector<VertexElement<2, 3>*> faces;
        unsigned num_faces_in_element = element_data.Faces.size();
        for (unsigned i = 0; i < num_faces_in_element; i++)
        {
            // Get the data for this face
            ElementData face_data = element_data.Faces[i];

            // Get the face index
            unsigned face_index = (unsigned)face_data.AttributeValue;

            // Get the nodes owned by this face
            std::vector<Node<3>*> nodes_in_face;
            unsigned num_nodes_in_face = face_data.NodeIndices.size();
            for (unsigned j = 0; j < num_nodes_in_face; j++)
            {
                assert(face_data.NodeIndices[j] < this->mNodes.size());
                nodes_in_face.push_back(this->mNodes[face_data.NodeIndices[j]]);
            }

            // If this face index is not already encountered, create a new face and update faces_counted...
            if (faces_counted.find(face_index) == faces_counted.end())
            {
                // Use nodes and index to construct this face
                VertexElement<2, 3>* p_face = new VertexElement<2, 3>(face_index, nodes_in_face);
                mFaces.push_back(p_face);
                faces_counted.insert(face_index);
                faces.push_back(p_face);
            }
            else
            {
                // ... otherwise use the member of mFaces with this index
                bool face_added = false;
                for (unsigned k = 0; k < mFaces.size(); k++)
                {
                    if (mFaces[k]->GetIndex() == face_index)
                    {
                        faces.push_back(mFaces[k]);
                        face_added = true;
                        break;
                    }
                }
                UNUSED_OPT(face_added);
                assert(face_added == true);
            }
        }

        ///\todo Store orientations? (#1076/#1377)
        std::vector<bool> orientations = std::vector<bool>(num_faces_in_element, true);

        // Use faces and index to construct this element
        VertexElement<3, 3>* p_element = new VertexElement<3, 3>(elem_index, faces, orientations, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(
    const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector;
    if (mpDelaunayMesh)
    {
        vector = mpDelaunayMesh->GetVectorFromAtoB(rLocationA, rLocationB);
    }
    else
    {
        vector = AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(rLocationA, rLocationB);
    }
    return vector;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVolumeOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double element_volume = 0.0;
    if (SPACE_DIM == 2)
    {
        // We map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
        c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

        unsigned num_nodes = p_element->GetNumNodes();
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
            c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            element_volume += 0.5 * (this_x * next_y - next_x * this_y);

            pos_1 = pos_2;
        }
    }
    else
    {
        // Loop over faces and add up pyramid volumes
        c_vector<double, SPACE_DIM> pyramid_apex = p_element->GetNodeLocation(0);
        for (unsigned face_index = 0; face_index < p_element->GetNumFaces(); face_index++)
        {
            // Get pointer to face
            VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = p_element->GetFace(face_index);

            // Calculate the area of the face and get unit normal to this face
            c_vector<double, SPACE_DIM> unit_normal;
            double face_area = CalculateUnitNormalToFaceWithArea(p_face, unit_normal);

            // Calculate the perpendicular distance from the plane of the face to the chosen apex
            c_vector<double, SPACE_DIM> base_to_apex = GetVectorFromAtoB(p_face->GetNodeLocation(0), pyramid_apex);
            double perpendicular_distance = fabs(inner_prod(base_to_apex, unit_normal));

            // Use these to calculate the volume of the pyramid formed by the face and the point pyramid_apex
            element_volume += face_area * perpendicular_distance / 3;
        }
    }

    // We take the absolute value just in case the nodes were really oriented clockwise
    return fabs(element_volume);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceAreaOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double surface_area = 0.0;
    if (SPACE_DIM == 2)
    {
        unsigned num_nodes = p_element->GetNumNodes();
        unsigned this_node_index = p_element->GetNodeGlobalIndex(0);
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            unsigned next_node_index = p_element->GetNodeGlobalIndex((local_index + 1) % num_nodes);

            surface_area += this->GetDistanceBetweenNodes(this_node_index, next_node_index);
            this_node_index = next_node_index;
        }
    }
    else
    {
        // Loop over faces and add up areas
        for (unsigned face_index = 0; face_index < p_element->GetNumFaces(); face_index++)
        {
            surface_area += CalculateAreaOfFace(p_element->GetFace(face_index));
        }
    }
    return surface_area;
}

//////////////////////////////////////////////////////////////////////
//                        2D-specific methods                       //
//////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time
    assert(ELEMENT_DIM == SPACE_DIM); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get the element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    // Initialise boolean
    bool element_includes_point = false;

    ///\todo (see #2387 and #2401) Investigate why the commented implementation causes Test2DVertexBasedCryptRepresentativeSimulation to fail
    //    // Initialise boolean
    //    bool element_includes_point = true;
    //
    //    unsigned winding_number = 0;
    //
    //    c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
    //    c_vector<double, SPACE_DIM> test_point = this->GetVectorFromAtoB(first_node_location, rTestPoint);
    //    c_vector<double, SPACE_DIM> this_node_location = zero_vector<double>(SPACE_DIM);
    //
    //    // Loop over edges of the element
    //    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    //    {
    //        c_vector<double, SPACE_DIM> untransformed_vector = p_element->GetNodeLocation((local_index+1)%num_nodes);
    //        c_vector<double, SPACE_DIM> next_node_location = this->GetVectorFromAtoB(first_node_location, untransformed_vector);
    //
    //        // If this edge is crossing upward...
    //        if (this_node_location[1] <= test_point[1])
    //        {
    //            if (next_node_location[1] > test_point[1])
    //            {
    //                double is_left =  (next_node_location[0] - this_node_location[0])*(test_point[1] - this_node_location[1])
    //                                 - (test_point[0] - this_node_location[0])*(next_node_location[1] - this_node_location[1]);
    //
    //                // ...and the test point is to the left of the edge...
    //                if (is_left > DBL_EPSILON)
    //                {
    //                    // ...then there is a valid upward edge-ray intersection to the right of the test point
    //                    winding_number++;
    //                }
    //            }
    //        }
    //        else
    //        {
    //            // ...otherwise, if the edge is crossing downward
    //            if (next_node_location[1] <= test_point[1])
    //            {
    //                double is_left =  (next_node_location[0] - this_node_location[0])*(test_point[1] - this_node_location[1])
    //                                 - (test_point[0] - this_node_location[0])*(next_node_location[1] - this_node_location[1]);
    //
    //                // ...and the test point is to the right of the edge...
    //                if (is_left < -DBL_EPSILON)
    //                {
    //                    // ...then there is a valid downward edge-ray intersection to the right of the test point
    //                    winding_number--;
    //                }
    //            }
    //        }
    //        this_node_location = next_node_location;
    //    }
    //
    //    if (winding_number == 0)
    //    {
    //        element_includes_point = false;
    //    }
    /////////////////////////////////////////////////////////

    // Remap the origin to the first vertex to allow alternative distance metrics to be used in subclasses
    c_vector<double, SPACE_DIM> first_vertex = p_element->GetNodeLocation(0);
    c_vector<double, SPACE_DIM> test_point = GetVectorFromAtoB(first_vertex, rTestPoint);

    // Loop over edges of the element
    c_vector<double, SPACE_DIM> vertexA = zero_vector<double>(SPACE_DIM);
    for (unsigned local_index = 0; local_index < num_nodes; local_index++)
    {
        // Check if this edge crosses the ray running out horizontally (increasing x, fixed y) from the test point
        c_vector<double, SPACE_DIM> vector_a_to_point = GetVectorFromAtoB(vertexA, test_point);

        // Pathological case - test point coincides with vertexA
        // (we will check vertexB next time we go through the for loop)
        if (norm_2(vector_a_to_point) < DBL_EPSILON)
        {
            return false;
        }

        c_vector<double, SPACE_DIM> vertexB = GetVectorFromAtoB(first_vertex, p_element->GetNodeLocation((local_index + 1) % num_nodes));
        c_vector<double, SPACE_DIM> vector_b_to_point = GetVectorFromAtoB(vertexB, test_point);
        c_vector<double, SPACE_DIM> vector_a_to_b = GetVectorFromAtoB(vertexA, vertexB);

        // Pathological case - ray coincides with horizontal edge
        if ((fabs(vector_a_to_b[1]) < DBL_EPSILON) && (fabs(vector_a_to_point[1]) < DBL_EPSILON) && (fabs(vector_b_to_point[1]) < DBL_EPSILON))
        {
            if ((vector_a_to_point[0] > 0) != (vector_b_to_point[0] > 0))
            {
                return false;
            }
        }

        // Non-pathological case
        // A and B on different sides of the line y = test_point[1]
        if ((vertexA[1] > test_point[1]) != (vertexB[1] > test_point[1]))
        {
            // Intersection of y=test_point[1] and vector_a_to_b is on the right of test_point
            if (test_point[0] < vertexA[0] + vector_a_to_b[0] * vector_a_to_point[1] / vector_a_to_b[1])
            {
                element_includes_point = !element_includes_point;
            }
        }

        vertexA = vertexB;
    }
    return element_includes_point;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time
    assert(ELEMENT_DIM == SPACE_DIM); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get the element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    double min_squared_normal_distance = DBL_MAX;
    unsigned min_distance_edge_index = UINT_MAX;

    // Loop over edges of the element
    for (unsigned local_index = 0; local_index < num_nodes; local_index++)
    {
        // Get the end points of this edge
        c_vector<double, SPACE_DIM> vertexA = p_element->GetNodeLocation(local_index);
        c_vector<double, SPACE_DIM> vertexB = p_element->GetNodeLocation((local_index + 1) % num_nodes);

        c_vector<double, SPACE_DIM> vector_a_to_point = this->GetVectorFromAtoB(vertexA, rTestPoint);
        c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);
        double distance_a_to_b = norm_2(vector_a_to_b);

        c_vector<double, SPACE_DIM> edge_ab_unit_vector = vector_a_to_b / norm_2(vector_a_to_b);
        double distance_parallel_to_edge = inner_prod(vector_a_to_point, edge_ab_unit_vector);

        double squared_distance_normal_to_edge = SmallPow(norm_2(vector_a_to_point), 2) - SmallPow(distance_parallel_to_edge, 2);

        /*
         * If the point lies almost bang on the supporting line of the edge, then snap to the line.
         * This allows us to do floating point tie-breaks when line is exactly at a node.
         * We adopt a similar approach if the point is at the same position as a point in the
         * element.
         */
        if (squared_distance_normal_to_edge < DBL_EPSILON)
        {
            squared_distance_normal_to_edge = 0.0;
        }

        if (fabs(distance_parallel_to_edge) < DBL_EPSILON)
        {
            distance_parallel_to_edge = 0.0;
        }
        else if (fabs(distance_parallel_to_edge - distance_a_to_b) < DBL_EPSILON)
        {
            distance_parallel_to_edge = distance_a_to_b;
        }

        // Make sure node is within the confines of the edge and is the nearest edge to the node \this breaks for convex elements
        if (squared_distance_normal_to_edge < min_squared_normal_distance && distance_parallel_to_edge >= 0 && distance_parallel_to_edge <= distance_a_to_b)
        {
            min_squared_normal_distance = squared_distance_normal_to_edge;
            min_distance_edge_index = local_index;
        }
    }

    assert(min_distance_edge_index < num_nodes);
    return min_distance_edge_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> VertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMomentsOfElement(unsigned index)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    // Define helper variables
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    unsigned num_nodes = p_element->GetNumNodes();
    c_vector<double, 3> moments = zero_vector<double>(3);

    // Since we compute I_xx, I_yy and I_xy about the centroid, we must shift each vertex accordingly
    c_vector<double, SPACE_DIM> centroid = GetCentroidOfElement(index);

    c_vector<double, SPACE_DIM> this_node_location = p_element->GetNodeLocation(0);
    c_vector<double, SPACE_DIM> pos_1 = this->GetVectorFromAtoB(centroid, this_node_location);

    for (unsigned local_index = 0; local_index < num_nodes; local_index++)
    {
        unsigned next_index = (local_index + 1) % num_nodes;
        c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation(next_index);
        c_vector<double, SPACE_DIM> pos_2 = this->GetVectorFromAtoB(centroid, next_node_location);

        double signed_area_term = pos_1(0) * pos_2(1) - pos_2(0) * pos_1(1);
        // Ixx
        moments(0) += (pos_1(1) * pos_1(1) + pos_1(1) * pos_2(1) + pos_2(1) * pos_2(1)) * signed_area_term;

        // Iyy
        moments(1) += (pos_1(0) * pos_1(0) + pos_1(0) * pos_2(0) + pos_2(0) * pos_2(0)) * signed_area_term;

        // Ixy
        moments(2) += (pos_1(0) * pos_2(1) + 2 * pos_1(0) * pos_1(1) + 2 * pos_2(0) * pos_2(1) + pos_2(0) * pos_1(1)) * signed_area_term;

        pos_1 = pos_2;
    }

    moments(0) /= 12;
    moments(1) /= 12;
    moments(2) /= 24;

    /*
     * If the nodes owned by the element were supplied in a clockwise rather
     * than anticlockwise manner, or if this arose as a result of enforcing
     * periodicity, then our computed quantities will be the wrong sign, so
     * we need to fix this.
     */
    if (moments(0) < 0.0)
    {
        moments(0) = -moments(0);
        moments(1) = -moments(1);
        moments(2) = -moments(2);
    }
    return moments;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetShortAxisOfElement(unsigned index)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);

    // Calculate the moments of the element about its centroid (recall that I_xx and I_yy must be non-negative)
    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    // Normalise the moments vector to remove problem of a very small discriminant (see #2874)
    moments /= norm_2(moments);

    // If the principal moments are equal...
    double discriminant = (moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2);
    if (fabs(discriminant) < DBL_EPSILON)
    {
        // ...then every axis through the centroid is a principal axis, so return a random unit vector
        short_axis(0) = RandomNumberGenerator::Instance()->ranf();
        short_axis(1) = sqrt(1.0 - short_axis(0) * short_axis(0));
    }
    else
    {
        // If the product of inertia is zero, then the coordinate axes are the principal axes
        if (fabs(moments(2)) < DBL_EPSILON)
        {
            if (moments(0) < moments(1))
            {
                short_axis(0) = 0.0;
                short_axis(1) = 1.0;
            }
            else
            {
                short_axis(0) = 1.0;
                short_axis(1) = 0.0;
            }
        }
        else
        {
            // Otherwise we find the eigenvector of the inertia matrix corresponding to the largest eigenvalue
            double lambda = 0.5 * (moments(0) + moments(1) + sqrt(discriminant));

            short_axis(0) = 1.0;
            short_axis(1) = (moments(0) - lambda) / moments(2);

            // Normalise the short axis before returning it
            short_axis /= norm_2(short_axis);
        }
    }

    return short_axis;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    unsigned num_nodes_in_element = pElement->GetNumNodes();
    unsigned next_local_index = (localIndex + 1) % num_nodes_in_element;

    // We add an extra num_nodes_in_element in the line below as otherwise this term can be negative, which breaks the % operator
    unsigned previous_local_index = (num_nodes_in_element + localIndex - 1) % num_nodes_in_element;

    c_vector<double, SPACE_DIM> previous_node_location = pElement->GetNodeLocation(previous_local_index);
    c_vector<double, SPACE_DIM> next_node_location = pElement->GetNodeLocation(next_local_index);
    c_vector<double, SPACE_DIM> difference_vector = this->GetVectorFromAtoB(previous_node_location, next_node_location);

    c_vector<double, SPACE_DIM> area_gradient;

    area_gradient[0] = 0.5 * difference_vector[1];
    area_gradient[1] = -0.5 * difference_vector[0];

    return area_gradient;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    unsigned num_nodes_in_element = pElement->GetNumNodes();

    // We add an extra num_nodes_in_element-1 in the line below as otherwise this term can be negative, which breaks the % operator
    unsigned previous_local_index = (num_nodes_in_element + localIndex - 1) % num_nodes_in_element;

    unsigned this_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned previous_global_index = pElement->GetNodeGlobalIndex(previous_local_index);

    double previous_edge_length = this->GetDistanceBetweenNodes(this_global_index, previous_global_index);
    assert(previous_edge_length > DBL_EPSILON);

    c_vector<double, SPACE_DIM> previous_edge_gradient = this->GetVectorFromAtoB(pElement->GetNodeLocation(previous_local_index), pElement->GetNodeLocation(localIndex)) / previous_edge_length;

    return previous_edge_gradient;
}

/*
//added 8-6-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetDampingConstant(double DampingConstant)
{
    mDampingConstant = DampingConstant;
}
*/

/*
//added 8-6-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetDampingConstant()
{
    return mDampingConstant;
}
*/


//added 7-29-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetRestLength(double RestLength)
{
    mRestLength = RestLength;
}



//added 7-29-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetSpringConstant(double SpringConstant)
{
    mSpringConstant = SpringConstant;
}


//added 8-11-23
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetSpringConstantDurotaxis(double SubstrateNodes, double minStiffness, double maxStiffness)
{
	//std::cout<<"minstiffness"<<minStiffness<<std::endl;
	
	//std::cout<<"this->mNodes.size()\t"<<this->mNodes.size()<<std::endl;
	//exit(0);
	std::vector<double> vect;
	//vect.push_back(this->mNodes.size());	
	double diff = (maxStiffness - minStiffness)/(this->mNodes.size() - 1);
	for (unsigned node_index = 0; node_index < this->mNodes.size(); node_index++)
    	{
		double k_s = minStiffness + (node_index*diff);
		//std::cout<<node_index+1-uv<<std::endl;
		Node<SPACE_DIM>* p_node = this->GetNode(node_index);
		p_node->AddElementStiffnessContribution(k_s);
		vect.push_back(k_s);
        	//std::cout<<node_index<<"\t"<<p_node->rGetElementStiffnessContribution()<<std::endl;
    	}
	//exit(0);
	
	for (unsigned element_index = 0; element_index < this->mElements.size(); element_index++)
    	{
		VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(element_index);
		double total_stiffness = 0;
		double mean_stiffness = 0;
		for (unsigned local_index = 0; local_index < p_element->GetNumNodes(); local_index++)
   		{
        		// Get a pointer to this node
        		Node<SPACE_DIM>* p_node = p_element->GetNode(local_index);
			total_stiffness = total_stiffness + p_node->rGetElementStiffnessContribution();
		
		}
		mean_stiffness = (total_stiffness/2);
		p_element->SetAttribute(mean_stiffness);
		//std::cout<<p_element->GetAttribute()<<std::endl;    	//for printing element stiffnesses, commented 8-13-23
	}
		
}



//added 12-22-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetVariableDampingConstant(double constantDampingConstant, double minDampingConstant, bool variableDampingConstant, double minNormalDistribution, double maxNormalDistribution, double stdevDampingConstant, double cellDampingFactor)
{
	double numcells = this->mNodes.size();
 	double uv = ((numcells-1)/2)+1;
	double s = stdevDampingConstant;
	//std::normal_distribution<double> distribution(a, 1.0);
	//std::vector<double> vect;
	//vect.push_back(this->mNodes.size());
	for (unsigned node_index = 0; node_index < this->mNodes.size(); node_index++)
    	{
		Node<SPACE_DIM>* p_node = this->GetNode(node_index);
		if (variableDampingConstant)
		{		
			//double d_c = minDampingConstant + 1*(cellDampingFactor)*(1/(s*sqrt(2*M_PI)))*exp(-(pow(((node_index+1-uv)/s),2.0))/2);

			double d_c = minDampingConstant + cellDampingFactor*(-minNormalDistribution + 1*(1)*(1/(s*sqrt(2*M_PI)))*exp(-(pow(((node_index+1-uv)/s),2.0))/2));
			//double d_c = cellDampingFactor*d_c;

			//std::cout<<node_index+1-uv<<std::endl;
			p_node->AddDampingConstant(d_c);
			//vect.push_back(k_s);
	        	std::cout<<"d_c\t"<<d_c<<std::endl;
		}
		else
		{
			p_node->AddDampingConstant(constantDampingConstant);
		}
	}
}

//added 11-25-21, modified 12-19-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetVariableSpringConstant(double CellNodes, double StdDev, double Factor)
{
    //mSpringConstant = SpringConstant;
	//double a = mElements.size()+1;
 	//VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
	//Node<SPACE_DIM>* p_node = p_element->GetNode(local_index);
	double uv = (CellNodes/2)+1;
	double s = StdDev;
	//std::normal_distribution<double> distribution(a, 1.0);
	std::vector<double> vect;
	//vect.push_back(this->mNodes.size());
	for (unsigned node_index = 0; node_index < this->mNodes.size(); node_index++)
    	{
		double k_s = 1*(1/(s*sqrt(2*M_PI)))*exp(-(pow(((node_index+1-uv)/s),2.0))/2);
		//std::cout<<node_index+1-uv<<std::endl;
		Node<SPACE_DIM>* p_node = this->GetNode(node_index);
		p_node->AddElementStiffnessContribution(k_s);
		vect.push_back(k_s);
        	//std::cout<<p_node->rGetElementStiffnessContribution()<<std::endl;
    	}
	double max = *max_element(vect.begin(), vect.end());
	double min = *min_element(vect.begin(), vect.end());	
	//std::cout<<*max_element(vect.begin(), vect.end())<<std::endl;

	for (unsigned node_index = 0; node_index < this->mNodes.size(); node_index++)
    	{
		//double k_s = 1*(1/(s*sqrt(2*M_PI)))*exp(-(pow(((node_index+1-uv)/s),2.0))/2);
		Node<SPACE_DIM>* p_node = this->GetNode(node_index);
		double stiffness = p_node->rGetElementStiffnessContribution();
		//stiffness = stiffness*(-1) + max + min + 5;
		stiffness = stiffness*(-1) + max + min;
		stiffness = Factor*(stiffness);
		p_node->AddElementStiffnessContribution(stiffness);
		//vect.push_back(k_s);
        	//std::cout<<p_node->rGetElementStiffnessContribution()<<std::endl;
    	}
	
	for (unsigned element_index = 0; element_index < this->mElements.size(); element_index++)
    	{
		VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(element_index);
		double total_stiffness = 0;
		double mean_stiffness = 0;
		for (unsigned local_index = 0; local_index < p_element->GetNumNodes(); local_index++)
   		{
        		// Get a pointer to this node
        		Node<SPACE_DIM>* p_node = p_element->GetNode(local_index);
			total_stiffness = total_stiffness + p_node->rGetElementStiffnessContribution();
			
		}
		mean_stiffness = (total_stiffness/2)*100;
		p_element->SetAttribute(mean_stiffness);
		std::cout<<p_element->GetAttribute()<<std::endl;
	}
		
	
}



//added 11-29-20
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNextEdgeVariableHookeSpringForceAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex)
{
	//double rest_length = 1/sqrt(3);

	//double rest_length = 0.4;

	double rest_length = mRestLength;	//added 7-29-21

    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

	//std::cout<<"before_next_local_index"<<std::endl;

    unsigned next_local_index = (localIndex + 1) % (pElement->GetNumNodes());

	//std::cout<<"after_next_local_index"<<std::endl;

    unsigned this_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned next_global_index = pElement->GetNodeGlobalIndex(next_local_index);

	//std::cout<<"after_next_global_index\t"<<next_global_index<<std::endl;


	//double spring_constant = mSpringConstant;   //modified 7-29-21

	double spring_constant = pElement->GetAttribute();   //added 11-29-21

    double next_edge_length = this->GetDistanceBetweenNodes(this_global_index, next_global_index);
    assert(next_edge_length > DBL_EPSILON);
	double overlap = next_edge_length - rest_length;
		

    c_vector<double, SPACE_DIM> next_edge_gradient = spring_constant*overlap*this->GetVectorFromAtoB(pElement->GetNodeLocation(next_local_index), pElement->GetNodeLocation(localIndex)) / next_edge_length;


	//std::cout<<"element spring constant\t"<<pElement->GetAttribute()<<"\t"<<"overlap\t"<<overlap<<"\t"<<"next_edge_gradient\t"<<next_edge_gradient<<std::endl;	//added 8-13-23

	//std::cout<<"this->GetVectorFromAtoB(pElement->GetNodeLocation(next_local_index), pElement->GetNodeLocation(localIndex)) / next_edge_length\t"<<this->GetVectorFromAtoB(pElement->GetNodeLocation(next_local_index), pElement->GetNodeLocation(localIndex)) / next_edge_length<<std::endl;	//added 10/10/22 for figuring out units, also this measures whether force is negative or positive, the magnitude is always 1
	
	//std::cout<<"overlap@hookeforce\t"<<overlap<<"\tlocalindex\t"<<localIndex<<"\tnext_edge_grad\t"<<next_edge_gradient<<std::endl;	

    return next_edge_gradient;
}




//added 8-6-20	modified 12-3-21
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNextEdgeHookeSpringForceAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex)
{
	//double rest_length = 1/sqrt(3);

	//double rest_length = 0.4;

	double rest_length = mRestLength;	//added 7-29-21

    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

	//std::cout<<"before_next_local_index"<<std::endl;

    unsigned next_local_index = (localIndex + 1) % (pElement->GetNumNodes());

	//std::cout<<"after_next_local_index"<<std::endl;

    unsigned this_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned next_global_index = pElement->GetNodeGlobalIndex(next_local_index);

	//std::cout<<"after_next_global_index\t"<<next_global_index<<std::endl;


	double spring_constant = mSpringConstant;   //modified 7-29-21	

    double next_edge_length = this->GetDistanceBetweenNodes(this_global_index, next_global_index);
    assert(next_edge_length > DBL_EPSILON);
	double overlap = next_edge_length - rest_length;

	if (abs(overlap) < 0.000000001)		//added 12-8-21
	{
		overlap = 0;
        }
//	std::cout<<"next_edge_length\t"<<next_edge_length<<"\trest_length\t"<<rest_length<<"\toverlap\t"<<overlap<<std::endl;	

    c_vector<double, SPACE_DIM> next_edge_gradient = spring_constant*overlap*this->GetVectorFromAtoB(pElement->GetNodeLocation(next_local_index), pElement->GetNodeLocation(localIndex)) / next_edge_length;

	//std::cout<<"overlap@hookeforce\t"<<overlap<<"\tlocalindex\t"<<localIndex<<"\tnext_edge_grad\t"<<next_edge_gradient<<std::endl;	

    return next_edge_gradient;
}





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    unsigned next_local_index = (localIndex + 1) % (pElement->GetNumNodes());

    unsigned this_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned next_global_index = pElement->GetNodeGlobalIndex(next_local_index);

    double next_edge_length = this->GetDistanceBetweenNodes(this_global_index, next_global_index);
    assert(next_edge_length > DBL_EPSILON);

    c_vector<double, SPACE_DIM> next_edge_gradient = this->GetVectorFromAtoB(pElement->GetNodeLocation(next_local_index), pElement->GetNodeLocation(localIndex)) / next_edge_length;

    return next_edge_gradient;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE

    c_vector<double, SPACE_DIM> previous_edge_gradient = GetPreviousEdgeGradientOfElementAtNode(pElement, localIndex);
    c_vector<double, SPACE_DIM> next_edge_gradient = GetNextEdgeGradientOfElementAtNode(pElement, localIndex);

    return previous_edge_gradient + next_edge_gradient;
}

//////////////////////////////////////////////////////////////////////
//                        3D-specific methods                       //
//////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateUnitNormalToFaceWithArea(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, c_vector<double, SPACE_DIM>& rNormal)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // As we are in 3D, the face must have at least three vertices
    assert(pFace->GetNumNodes() >= 3u);

    // Reset the answer
    rNormal = zero_vector<double>(SPACE_DIM);

    c_vector<double, SPACE_DIM> v_minus_v0 = this->GetVectorFromAtoB(pFace->GetNode(0)->rGetLocation(), pFace->GetNode(1)->rGetLocation());
    for (unsigned local_index = 2; local_index < pFace->GetNumNodes(); local_index++)
    {
        c_vector<double, SPACE_DIM> vnext_minus_v0 = this->GetVectorFromAtoB(pFace->GetNode(0)->rGetLocation(), pFace->GetNode(local_index)->rGetLocation());
        rNormal += VectorProduct(v_minus_v0, vnext_minus_v0);
        v_minus_v0 = vnext_minus_v0;
    }
    double magnitude = norm_2(rNormal);
    if (magnitude != 0.0)
    {
        // Normalize the normal vector
        rNormal /= magnitude;
        // If all points are co-located, then magnitude==0.0 and there is potential for a floating point exception
        // here if we divide by zero, so we'll move on.
    }
    return magnitude / 2.0;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateAreaOfFace(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get the unit normal to the plane of this face
    c_vector<double, SPACE_DIM> unit_normal;
    return CalculateUnitNormalToFaceWithArea(pFace, unit_normal);
}

/// Specialization to avoid compiler error about zero-sized arrays
#if defined(__xlC__)
template <>
double VertexMesh<1, 1>::CalculateAreaOfFace(VertexElement<0, 1>* pFace)
{
    NEVER_REACHED;
}
#endif

// Explicit instantiation
template class VertexMesh<1, 1>;
template class VertexMesh<1, 2>;
template class VertexMesh<1, 3>;
template class VertexMesh<2, 2>;
template class VertexMesh<2, 3>;
template class VertexMesh<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh)
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

#include "SubstrateHookeForce_ab.hpp"

template<unsigned DIM>
SubstrateHookeForce_ab<DIM>::SubstrateHookeForce_ab()
   : AbstractForce<DIM>(),
     mNagaiHondaDeformationEnergyParameter(100.0), // This is 1.0 in the Nagai & Honda paper.
     mNagaiHondaMembraneSurfaceEnergyParameter(10.0), // This is 0.1 in the Nagai & Honda paper.
     mNagaiHondaCellCellAdhesionEnergyParameter(0.5), // This corresponds to a value of 1.0 for
                                                      // the sigma parameter in the Nagai & Honda
                                                      // paper. In the paper, the sigma value is
                                                      // set to 0.01.
     mNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0), // This is 0.01 in the Nagai & Honda paper.
	mSubstrateHookeForceCellCellAdhesionEnergyParameter(0.5), //added 8-5-20
	mSubstrateHookeForceCellBoundaryAdhesionEnergyParameter(1.0) //added 8-5-20
{
}

template<unsigned DIM>
SubstrateHookeForce_ab<DIM>::~SubstrateHookeForce_ab()
{
}

//added 8-4-20

template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::AddForceContribution_ab1(MutableVertexMesh<DIM>& rMesh_ab)
{
    //std::cout<<"inside SubstrateHookeForce_ab"<<std::endl;
    //exit(0);
    // Throw an exception message if not using a VertexBasedCellPopulation
	MutableVertexMesh<DIM>* p_mesh_ab = static_cast<MutableVertexMesh<DIM>*>(&rMesh_ab);
	unsigned num_nodes = p_mesh_ab->GetNumNodes();
	//std::cout<<"p_mesh_ab->GetNumNodes()\t"<<p_mesh_ab->GetNumNodes()<<std::endl;

	for (unsigned node_index=0; node_index<num_nodes; node_index++)
   	{
	        
		Node<DIM>* p_this_node = p_mesh_ab->GetNode(node_index);
		//Node<DIM>* p_this_node = p_mesh_ab->GetNode(3);
		// Find the indices of the elements owned by this node

		c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);

        	std::set<unsigned> containing_elem_indices = p_this_node->rGetContainingElementIndices();
		//std::cout<<"std::set<unsigned> containing_elem_indices\t"<<containing_elem_indices.size()<<std::endl;
		//std::cout<<"std::set<unsigned> containing_elem_indices\t"<<containing_elem_indices[0]<<"\t"<<containing_elem_indices[1]<<std::endl;
		// Iterate over these elements
		//unsigned abc = 1;
        	for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             	iter != containing_elem_indices.end();
             	++iter)
        	{
			//std::cout<<"abc\t"<<abc<<std::endl;			
			// Get this element, its index and its number of nodes
            		VertexElement<DIM, DIM>* p_element = p_mesh_ab->GetElement(*iter);
            		//unsigned elem_index = p_element->GetIndex();
			//std::cout<<"elem_index\t"<<elem_index<<std::endl;
            		unsigned num_nodes_elem = p_element->GetNumNodes();
			//std::cout<<"num_nodes_elem\t"<<num_nodes_elem<<std::endl;
			// Find the local index of this node in this element
            		unsigned local_index = p_element->GetNodeLocalIndex(node_index);
			
			//std::cout<<"element_stiffness\t"<<p_element->GetAttribute()<<std::endl;		//added 8-13-23 to check if durotaxis stiffnesses were passed corrently or not.
			//std::cout<<"local_index\t"<<local_index<<std::endl;
			//exit(0);
			// Get the previous and next nodes in this element
           		unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            		Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
			//std::cout<<"previous_node_local_index\t"<<previous_node_local_index<<std::endl;
			unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            		Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);
			//std::cout<<"next_node_local_index\t"<<next_node_local_index<<std::endl;
			//abc++;
			// Compute the adhesion parameter for each of these edges
            		double previous_edge_adhesion_parameter = GetAdhesionParameter_ab(p_previous_node, p_this_node, *p_mesh_ab);
            		double next_edge_adhesion_parameter = GetAdhesionParameter_ab(p_this_node, p_next_node, *p_mesh_ab);
			//std::cout<<"abcdgejjg"<<std::endl;
		//	std::cout<<"previous_edge_adhesion_parameter\t"<<previous_edge_adhesion_parameter<<"\t"<<"next_edge_adhesion_parameter\t"<<next_edge_adhesion_parameter<<std::endl;
			//std::cout<<"sub_node_index_is\t"<<node_index<<"\tlocal_index_is\t"<<local_index<<std::endl;

			//double spring_constant = mSpringConstant;	//added 12-3-21 

			//double rest_length = 

			//c_vector<double, DIM> previous_edge_force = -(p_mesh_ab->GetNextEdgeHookeSpringForceAtNode(p_element, previous_node_local_index));	//modified 12-3-21	//commented 8-13-23 for implementing durotaxis
			
			c_vector<double, DIM> previous_edge_force = -(p_mesh_ab->GetNextEdgeVariableHookeSpringForceAtNode(p_element, previous_node_local_index));	//added 8-13-23 for durotaxis

			//c_vector<double, DIM> next_edge_force = p_mesh_ab->GetNextEdgeHookeSpringForceAtNode(p_element, local_index);	//modified 12-3-21	commented 8-13-23 for implementing durotaxis
			
			c_vector<double, DIM> next_edge_force = p_mesh_ab->GetNextEdgeVariableHookeSpringForceAtNode(p_element, local_index);		//added 8-13-23 for durotaxis

			//std::cout<<"sub_node_index_is\t"<<node_index<<std::endl;

			adhesion_contribution -= previous_edge_adhesion_parameter*previous_edge_force + next_edge_adhesion_parameter*next_edge_force;
			//std::cout<<"abcdgejjg8888888888"<<std::endl;
		}
		//std::cout<<"insidehooke"<<std::endl;
		//std::cout<<"sub_node_index_is\t"<<node_index<<"\tforce_b4_hooke_is\t"<<(p_mesh_ab->GetNode(node_index)->rGetAppliedForce())[0]<<"\t"<<(p_mesh_ab->GetNode(node_index)->rGetAppliedForce())[1]<<std::endl;

	//	std::cout<<"substratehookeforce_on_node[0]\t"<<adhesion_contribution[0]<<"\t"<<"substratehookeforce_on_node[1]\t"<<adhesion_contribution[1]<<std::endl;

		p_mesh_ab->GetNode(node_index)->AddAppliedForceContribution(adhesion_contribution);

		//std::cout<<"sub_node_index_is\t"<<node_index<<"\tforce_after_hooke_is\t"<<(p_mesh_ab->GetNode(node_index)->rGetAppliedForce())[0]<<"\t"<<(p_mesh_ab->GetNode(node_index)->rGetAppliedForce())[1]<<"\n"<<std::endl;
	}
//exit(0);
}

//added 8-5-20
template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetAdhesionParameter_ab(Node<DIM>* pNodeA, Node<DIM>* pNodeB, MutableVertexMesh<DIM>& rMesh_ab)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    double adhesion_parameter = GetSubstrateHookeForceCellCellAdhesionEnergyParameter();

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = GetSubstrateHookeForceCellBoundaryAdhesionEnergyParameter();
    }

    return adhesion_parameter;
}

//added 8-5-20
template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetSubstrateHookeForceCellCellAdhesionEnergyParameter()
{
    return mSubstrateHookeForceCellCellAdhesionEnergyParameter;
}

//added 8-5-20
template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetSubstrateHookeForceCellBoundaryAdhesionEnergyParameter()
{
    return mSubstrateHookeForceCellBoundaryAdhesionEnergyParameter;
}


//added 6-25-20
template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::AddForceContribution_ab(AbstractCellPopulation<DIM>& rCellPopulation, MutableVertexMesh<DIM>& rMesh_ab)
{
    std::cout<<"tractor"<<std::endl;
 //added 6-25-20
/*   for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = rMesh_ab.GetNodeIteratorBegin();
             node_iter != rMesh_ab.GetNodeIteratorEnd();
             ++node_iter)
    {
             std::cout<<"MIAsma"<<std::endl;            
             //old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
    }*/

	//}
    //out<<"loveandhate676767"<<std::endl;
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SubstrateHookeForce_ab is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    unsigned num_nodes = p_cell_population->GetNumNodes();
std::cout<<"tractor22"<<std::endl;

//std::cout<<"tractor33"<<std::endl;
    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
	     std::cout<<"tractor44"<<std::endl;   
	Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
		//std::cout<<"tractor55"<<std::endl;
	if (p_this_node->IsFocal()==true)
	{	//std::cout<<"tractor66"<<std::endl;	
		std::cout<<"p_this_node->rGetSubstrateNodeIndex()\t"<<p_this_node->rGetSubstrateNodeIndex()<<std::endl;									
		Node<DIM>* p_sub_node = rMesh_ab.GetNode(p_this_node->rGetSubstrateNodeIndex());	
		const c_vector<double, DIM>& r_cell_location = p_this_node->rGetLocation();
   		const c_vector<double, DIM>& r_sub_location = p_sub_node->rGetLocation();
		c_vector<double, DIM> unit_difference = GetVectorFromAtoB_ab(r_cell_location, r_sub_location);
		double distance_between_nodes = norm_2(unit_difference);
		//unit_difference /= distance_between_nodes;
		//std::cout<<"pakoda"<<std::endl;
		//std::cout<<unit_difference[0]<<distance_between_nodes<<std::endl;
		double rest_length = p_this_node->rGetSubstrateNodeDistance();
		double overlap = distance_between_nodes - rest_length;
		double spring_stiffness = 100;
		//std::cout<<"pakoda$$$"<<std::endl;
		//std::cout<<overlap<<std::endl;
		
		c_vector<double, DIM> force_on_node1 = spring_stiffness * unit_difference * overlap;

		std::cout<<"substratehookeforce_on_node[0]\t"<<force_on_node1[0]<<"\t"<<"substratehookeforce_on_node[1]\t"<<force_on_node1[1]<<std::endl;

		p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node1);
		rMesh_ab.GetNode(p_this_node->rGetSubstrateNodeIndex())->AddAppliedForceContribution(-force_on_node1);

	}
 
    }
}

template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    double adhesion_parameter = GetNagaiHondaCellCellAdhesionEnergyParameter();

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
    }

    return adhesion_parameter;
}

template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetNagaiHondaDeformationEnergyParameter()
{
    return mNagaiHondaDeformationEnergyParameter;
}

//added 6-30-20
template <unsigned DIM>
c_vector<double, DIM> SubstrateHookeForce_ab<DIM>::GetVectorFromAtoB_ab(
    const c_vector<double, DIM>& rLocationA, const c_vector<double, DIM>& rLocationB)
{
    c_vector<double, DIM> vector = rLocationB - rLocationA;
    return vector;
}

//added 6-30-20
template <unsigned DIM>
unsigned SubstrateHookeForce_ab<DIM>::GetNearestSubstrateNode2CellNode(const c_vector<double, DIM>& rLocationA, MutableVertexMesh<DIM>& rrMesh_ab)
{
    	// Hold the best distance from node to point found so far
    	// and the (local) node at which this was recorded
	unsigned best_node_index = 0u;
	double best_node_point_distance = DBL_MAX;

    	
   	// Now loop through the nodes, calculating the distance and updating best_node_point_distance
	for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = rrMesh_ab.GetNodeIteratorBegin();
             node_iter != rrMesh_ab.GetNodeIteratorEnd();
             ++node_iter)
	{
             	double node_point_distance = norm_2(GetVectorFromAtoB_ab((node_iter)->rGetLocation(), rLocationA));
		if (node_point_distance < best_node_point_distance)
	        {
		            best_node_index = (node_iter)->GetIndex();
		            best_node_point_distance = node_point_distance;
	        }	     	     
			
		//std::cout<<"node_point_distance"<<std::endl;	     
		//std::cout<<node_point_distance<<std::endl;            
        }
	//std::cout<<"best_node_point_distance"<<std::endl;	     
	//std::cout<<best_node_point_distance<<std::endl;   
	return best_node_index;

}



template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetNagaiHondaCellCellAdhesionEnergyParameter()
{
    return mNagaiHondaCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double SubstrateHookeForce_ab<DIM>::GetNagaiHondaCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::SetNagaiHondaDeformationEnergyParameter(double deformationEnergyParameter)
{
    mNagaiHondaDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::SetNagaiHondaMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mNagaiHondaCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}


//moved 6-30-20



template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    //std::cout<<"loveandhate266220000000hhhhhhhhhhhhh"<<std::endl;
    //out<<"loveandhate676767"<<std::endl;
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SubstrateHookeForce_ab is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use SubstrateHookeForce_ab");
        }
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
	        
	Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of three
         * parts - a cell deformation energy, a membrane surface tension energy
         * and an adhesion energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            c_vector<double, DIM> element_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            deformation_contribution -= 2*GetNagaiHondaDeformationEnergyParameter()*(element_areas[elem_index] - target_areas[elem_index])*element_area_gradient;

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the adhesion parameter for each of these edges
            double previous_edge_adhesion_parameter = GetAdhesionParameter(p_previous_node, p_this_node, *p_cell_population);
            double next_edge_adhesion_parameter = GetAdhesionParameter(p_this_node, p_next_node, *p_cell_population);

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, DIM> previous_edge_gradient = -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary adhesion (note the minus sign)
            adhesion_contribution -= previous_edge_adhesion_parameter*previous_edge_gradient + next_edge_adhesion_parameter*next_edge_gradient;

            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            c_vector<double, DIM> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
            double cell_target_perimeter = 2*sqrt(M_PI*target_areas[elem_index]);
            membrane_surface_tension_contribution -= 2*GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;
        }

        c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution + adhesion_contribution;
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
    }
}

template<unsigned DIM>
void SubstrateHookeForce_ab<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellCellAdhesionEnergyParameter>" << mNagaiHondaCellCellAdhesionEnergyParameter << "</NagaiHondaCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaCellBoundaryAdhesionEnergyParameter << "</NagaiHondaCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SubstrateHookeForce_ab<1>;
template class SubstrateHookeForce_ab<2>;
template class SubstrateHookeForce_ab<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SubstrateHookeForce_ab)

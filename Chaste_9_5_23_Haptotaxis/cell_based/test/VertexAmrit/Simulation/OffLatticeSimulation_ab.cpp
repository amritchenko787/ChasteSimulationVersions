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

#include "OffLatticeSimulation_ab.hpp"

#include <boost/make_shared.hpp>

#include "CellBasedEventHandler.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "ForwardEulerNumericalMethod_ab.hpp"
#include "StepSizeException.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::OffLatticeSimulation_ab(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
						MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>& rMesh_ab,                                                
						bool deleteCellPopulationInDestructor,
                                                bool initialiseCells)
    : AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
      mrMesh_ab(rMesh_ab) 	
{
    if (!dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OffLatticeSimulation_abs require a subclass of AbstractOffLatticeCellPopulation.");
    }
/*
    //{ added 6-17-20
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mrMesh_ab.GetNodeIteratorBegin();
             node_iter != mrMesh_ab.GetNodeIteratorEnd();
             ++node_iter)
        {
		std::cout<<"He-Man"<<std::endl;            
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }
	//exit(0);
*/
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>     //added 6-11-20
MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>& OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::rGetMesh_ab()
{
    return mrMesh_ab;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>     //added 6-11-20
const MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>& OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::rGetMesh_ab() const
{
    return mrMesh_ab;
}

//added 8-3-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::AddForce_ab(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce_ab)
{
    mForceCollection_ab.push_back(pForce_ab);
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::AddForce(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::RemoveAllForces()
{
    mForceCollection.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellPopulationBoundaryConditions()
{
    mBoundaryConditions.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::SetNumericalMethod(boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > pNumericalMethod)
{
    mpNumericalMethod = pNumericalMethod;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::GetNumericalMethod() const
{
    return mpNumericalMethod;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >& OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::rGetForceCollection() const
{
    return mForceCollection;
}

//added 7-31-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM, SPACE_DIM>::CloseVtkMetaFile_ab()
{
#ifdef CHASTE_VTK
    *(this->mpVtkMetaFile_ab) << "    </Collection>\n";
    *(this->mpVtkMetaFile_ab) << "</VTKFile>\n";
    (this->mpVtkMetaFile_ab)->close();
#endif //CHASTE_VTK
}




//added 7-31-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM, SPACE_DIM>::InitiateVtkMetaFile_ab(OutputFileHandler& routput_file_handler_ab)
{
#ifdef CHASTE_VTK
    this->mpVtkMetaFile_ab = routput_file_handler_ab.OpenOutputFile("results_ab.pvd");
    *(this->mpVtkMetaFile_ab) << "<?xml version=\"1.0\"?>\n";
    *(this->mpVtkMetaFile_ab) << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *(this->mpVtkMetaFile_ab) << "    <Collection>\n";
#endif //CHASTE_VTK
}




//added 7-29-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM, SPACE_DIM>::WriteVtkResultsToFile_ab(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
	//std::cout<<"why?"<<std::endl;
    // Create mesh writer for VTK output
    VertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(rDirectory, "results", false);

    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    mesh_writer.WriteVtkUsingMesh(mrMesh_ab, time.str());

    *(this->mpVtkMetaFile_ab) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile_ab) << num_timesteps;
    *(this->mpVtkMetaFile_ab) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile_ab) << num_timesteps;
    *(this->mpVtkMetaFile_ab) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}




//added 7-14-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::UpdateAdhesions()
{
	//VertexBasedCellPopulation<ELEMENT_DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<ELEMENT_DIM>*>(&mrCellPopulation);
    	//unsigned num_nodes = p_cell_population->GetNumNodes();
	double m = 1000000;
	double kon = 10;
	double koff = 10;
	//std::cout<<"not-joking1"<<std::endl;
	std::cout<<"@offlaattice_updateadhesion"<<std::endl;
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
		
		double random = (RandomNumberGenerator::Instance()->randMod(99999))/m;
		//std::cout<<"not-joking2"<<std::endl;            
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
		unsigned cell_node_index = (node_iter)->GetIndex();
		//std::cout<<"cell_node_index\t"<<cell_node_index<<std::endl; 
		Node<SPACE_DIM>* p_cell_node = this->mrCellPopulation.rGetMesh().GetNode(cell_node_index);
		if (p_cell_node->IsFocal()==false && random < (this->mDt)*kon)
		{	
			c_vector<double, SPACE_DIM> node_location_cell = p_cell_node->rGetLocation();
			//std::cout<<"kon"<<std::endl;		
			//std::cout<<node_location_cell[0]<<std::endl;
			unsigned nearest_sub_node_index = GetNearestSubstrateNode2CellNode(node_location_cell,cell_node_index);
			c_vector<double, SPACE_DIM> nearest_node_location_sub = (mrMesh_ab.GetNode(nearest_sub_node_index))->rGetLocation();
			//std::cout<<this->mDt<<std::endl;
			
			p_cell_node->HasFocal(true);
			p_cell_node->AddSubstrateNodeLocation(nearest_node_location_sub);
			//bool abc = p_cell_node->IsFocal();
			//std::cout<<abc<<std::endl;
			p_cell_node->AddSubstrateNodeDistance((mrMesh_ab.GetNode(nearest_sub_node_index))->rGetSubstrateNodeDistance());
			p_cell_node->AddSubstrateNodeIndex(nearest_sub_node_index);
			//std::cout<<"@offlaattice_updateadhesion"<<std::endl;
			std::cout<<"cell node index is\t"<<cell_node_index<<"\t"<<"nearest_sub_node_index is"<<nearest_sub_node_index<<"\t"<<"rgetsubstratenodeindex gives\t"<<p_cell_node->rGetSubstrateNodeIndex()<<std::endl;
			//std::cout<<p_cell_node->rGetSubstrateNodeDistance()<<std::endl;
		}
		else if (p_cell_node->IsFocal()==true && random < (this->mDt)*koff)
		{
			//std::cout<<"koff"<<std::endl;
			p_cell_node->HasFocal(false);
			p_cell_node->ClearNodeLocation();
			//see if you need to remove 'substratenodedistance' and 'substratenodeindex'. In order to do that you would need to store those first in objects like std::vector. First check wehter rewritting variable is more efficient than clearing ojects and then rewriting objects
			(mrMesh_ab.GetNode(p_cell_node->rGetSubstrateNodeIndex()))->HasFocal(false);
			(mrMesh_ab.GetNode(p_cell_node->rGetSubstrateNodeIndex()))->ClearNodeLocation();
		}
			


        }


}



//added 7-14-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::GetVectorFromAtoB_ab(
    const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector = rLocationB - rLocationA;
    return vector;
}

//added 7-14-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::GetNearestSubstrateNode2CellNode(const c_vector<double, SPACE_DIM>& rLocationA, unsigned& rcell_node_index)
{
    	// Hold the best distance from node to point found so far
    	// and the (local) node at which this was recorded
	unsigned best_node_index = 0u;
	double best_node_point_distance = DBL_MAX;

    	
   	// Now loop through the nodes, calculating the distance and updating best_node_point_distance
	for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = this->mrMesh_ab.GetNodeIteratorBegin();
             node_iter != this->mrMesh_ab.GetNodeIteratorEnd();
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
	
	Node<SPACE_DIM>* p_this_node1 = this->mrMesh_ab.GetNode(best_node_index);
	p_this_node1->AddSubstrateNodeLocation(rLocationA);
	p_this_node1->AddSubstrateNodeDistance(best_node_point_distance);
	p_this_node1->HasFocal(true);
	p_this_node1->AddSubstrateNodeIndex(rcell_node_index);
	//std::cout<<"best_node_point_distance"<<std::endl;
	//std::cout<<"begam"<<std::endl;		     
	//std::cout<<best_node_point_distance<<std::endl;   
	return best_node_index;

}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::UpdateSubstrateTopology()
{
    //CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    //std::cout<<"substratetopology"<<std::endl;
    double time_advanced_so_far = 0;
    double target_time_step  = this->mDt;
    double present_time_step = this->mDt;

    while (time_advanced_so_far < target_time_step)
    {
        /*
	//{ added 6-16-20
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mrMesh_ab.GetNodeIteratorBegin();
             node_iter != mrMesh_ab.GetNodeIteratorEnd();
             ++node_iter)
        {
		std::cout<<"He-Man"<<std::endl;            
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

	//}
        */
	// Store the initial node positions (these may be needed when applying boundary conditions)
        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations_ab;
	//std::cout<<"He-Man1111"<<std::endl;	
	//std::cout<<this<<std::endl;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrMesh_ab.GetNodeIteratorBegin();
             node_iter != this->mrMesh_ab.GetNodeIteratorEnd();
             ++node_iter)
        {
            old_node_locations_ab[&(*node_iter)] = (node_iter)->rGetLocation();
        }
//std::cout<<"He-Man2"<<std::endl;
        // Try to update node positions according to the numerical method
        try
        {
            mpNumericalMethod_ab->UpdateAllNodePositions_ab(present_time_step);
		//std::cout<<"after_UpdateAllNodePositions_ab"<<std::endl;
            //ApplyBoundaries(old_node_locations);
//std::cout<<"He-Man3"<<std::endl;
            // Successful time step! Update time_advanced_so_far
            time_advanced_so_far += present_time_step;

            // If using adaptive timestep, then increase the present_time_step (by 1% for now)
            if (mpNumericalMethod_ab->HasAdaptiveTimestep())
            {
                ///\todo #2087 Make this a settable member variable
                double timestep_increase = 0.01;
                present_time_step = std::min((1+timestep_increase)*present_time_step, target_time_step - time_advanced_so_far);
            }

        }
        catch (StepSizeException& e)
        {
            // Detects if a node has travelled too far in a single time step
            if (mpNumericalMethod_ab->HasAdaptiveTimestep())
            {
                // If adaptivity is switched on, revert node locations and choose a suitably smaller time step
                RevertToOldLocations(old_node_locations_ab);
                present_time_step = std::min(e.GetSuggestedNewStep(), target_time_step - time_advanced_so_far);
            }
            else
            {
                // If adaptivity is switched off, terminate with an error
                EXCEPTION(e.what());
            }
        }
    }

    //CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}







template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::UpdateCellLocationsAndTopology()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    //std::cout<<"He-Man1"<<std::endl;
    double time_advanced_so_far = 0;
    double target_time_step  = this->mDt;
    double present_time_step = this->mDt;

    while (time_advanced_so_far < target_time_step)
    {
        /*
	//{ added 6-16-20
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mrMesh_ab.GetNodeIteratorBegin();
             node_iter != mrMesh_ab.GetNodeIteratorEnd();
             ++node_iter)
        {
		std::cout<<"He-Man"<<std::endl;            
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

	//}
        */
	// Store the initial node positions (these may be needed when applying boundary conditions)
        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations;
	//std::cout<<"He-Man1111"<<std::endl;	
	//std::cout<<this<<std::endl;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }
//std::cout<<"He-Man2"<<std::endl;
        // Try to update node positions according to the numerical method
        try
        {
            mpNumericalMethod->UpdateAllNodePositions(present_time_step);
            ApplyBoundaries(old_node_locations);
//std::cout<<"He-Man3"<<std::endl;
            // Successful time step! Update time_advanced_so_far
            time_advanced_so_far += present_time_step;

            // If using adaptive timestep, then increase the present_time_step (by 1% for now)
            if (mpNumericalMethod->HasAdaptiveTimestep())
            {
                ///\todo #2087 Make this a settable member variable
                double timestep_increase = 0.01;
                present_time_step = std::min((1+timestep_increase)*present_time_step, target_time_step - time_advanced_so_far);
            }

        }
        catch (StepSizeException& e)
        {
            // Detects if a node has travelled too far in a single time step
            if (mpNumericalMethod->HasAdaptiveTimestep())
            {
                // If adaptivity is switched on, revert node locations and choose a suitably smaller time step
                RevertToOldLocations(old_node_locations);
                present_time_step = std::min(e.GetSuggestedNewStep(), target_time_step - time_advanced_so_far);
            }
            else
            {
                // If adaptivity is switched off, terminate with an error
                EXCEPTION(e.what());
            }
        }
    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::RevertToOldLocations(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > oldNodeLoctions)
{
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {
        (node_iter)->rGetModifiableLocation() = oldNodeLoctions[&(*node_iter)];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::ApplyBoundaries(std::map<Node<SPACE_DIM>*,c_vector<double, SPACE_DIM> > oldNodeLoctions)
{
    // Apply any boundary conditions
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        (*bcs_iter)->ImposeBoundaryCondition(oldNodeLoctions);
    }

    // Verify that each boundary condition is now satisfied
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        if (!((*bcs_iter)->VerifyBoundaryCondition()))
        {
            EXCEPTION("The cell population boundary conditions are incompatible.");
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::WriteVisualizerSetupFile()
{
	//std::cout<<"WriteVisualizerSetupFile111"<<std::endl;
    if (PetscTools::AmMaster())
    {
        for (unsigned i=0; i<this->mForceCollection.size(); i++)
        {
            this->mForceCollection[i]->WriteDataToVisualizerSetupFile(this->mpVizSetupFile);
        }

        this->mrCellPopulation.WriteDataToVisualizerSetupFile(this->mpVizSetupFile);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::SetupSolve()
{
   /**/ // Clear all forces in substrate mesh (added 8-4-20)
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrMesh_ab.GetNodeIteratorBegin();
         node_iter != this->mrMesh_ab.GetNodeIteratorEnd();
         ++node_iter)
    {
	std::cout<<"insidesetupsolve"<<std::endl;
        node_iter->ClearAppliedForce();
    }

	//(added 8-4-20)
	// Use a forward Euler method by default, unless a numerical method has been specified already
    if (mpNumericalMethod_ab == nullptr)
    {
        mpNumericalMethod_ab = boost::make_shared<ForwardEulerNumericalMethod_ab<ELEMENT_DIM, SPACE_DIM> >();
    }
	mpNumericalMethod_ab->SetForceCollection_ab(&mForceCollection_ab);//added 8-7-20
	mpNumericalMethod_ab->SetSubstrateMesh(&mrMesh_ab);//added 8-7-20


    // Clear all forces
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

    // Use a forward Euler method by default, unless a numerical method has been specified already
    if (mpNumericalMethod == nullptr)
    {
        mpNumericalMethod = boost::make_shared<ForwardEulerNumericalMethod_ab<ELEMENT_DIM, SPACE_DIM> >();
    }
    
    mpNumericalMethod->SetCellPopulation(dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation)));
    mpNumericalMethod->SetForceCollection(&mForceCollection); 
    //mpNumericalMethod->SetForceCollection_ab(&mForceCollection, &mForceCollection_ab);//added 8-4-20
    mpNumericalMethod->SetSubstrateMesh(&mrMesh_ab);//added 6-21-20

    

//std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > > mForceCollection;
//void SetForceCollection(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces);

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over forces
    *rParamsFile << "\n\t<Forces>\n";
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        // Output force details
        (*iter)->OutputForceInfo(rParamsFile);
    }
    *rParamsFile << "\t</Forces>\n";

    // Loop over cell population boundary conditions
    *rParamsFile << "\n\t<CellPopulationBoundaryConditions>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mBoundaryConditions.begin();
         iter != mBoundaryConditions.end();
         ++iter)
    {
        // Output cell boundary condition details
        (*iter)->OutputCellPopulationBoundaryConditionInfo(rParamsFile);
    }
    *rParamsFile << "\t</CellPopulationBoundaryConditions>\n";

    // Output numerical method details
    *rParamsFile << "\n\t<NumericalMethod>\n";
    mpNumericalMethod->OutputNumericalMethodInfo(rParamsFile);
    *rParamsFile << "\t</NumericalMethod>\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation_ab<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(rParamsFile);
}

// Explicit instantiation
template class OffLatticeSimulation_ab<1,1>;
template class OffLatticeSimulation_ab<1,2>;
template class OffLatticeSimulation_ab<2,2>;
template class OffLatticeSimulation_ab<1,3>;
template class OffLatticeSimulation_ab<2,3>;
template class OffLatticeSimulation_ab<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulation_ab)

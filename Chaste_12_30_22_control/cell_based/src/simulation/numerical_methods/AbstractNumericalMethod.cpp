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

#include "AbstractNumericalMethod.hpp"
#include "StepSizeException.hpp"
#include "Warnings.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellBasedEventHandler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::AbstractNumericalMethod()
    : mpCellPopulation(nullptr),
      mpForceCollection(nullptr),
      mUseAdaptiveTimestep(false),
      mUseUpdateNodeLocation(false),
      mGhostNodeForcesEnabled(true)
{
    // mpCellPopulation and mpForceCollection are initialized by the OffLatticeSimulation constructor
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::~AbstractNumericalMethod()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetCellPopulation(AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>* pPopulation)
{
    mpCellPopulation = pPopulation;

    // Set other member variables according to the type of the cell population
    if (dynamic_cast<NodeBasedCellPopulationWithBuskeUpdate<SPACE_DIM>*>(mpCellPopulation))
    {
        mUseUpdateNodeLocation = true;
        WARNING("Non-Euler steppers are not yet implemented for NodeBasedCellPopulationWithBuskeUpdate");
    }

    if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(mpCellPopulation))
    {
        mGhostNodeForcesEnabled = true;
    }
    else
    {
        mGhostNodeForcesEnabled = false;
    }
}

//added 7-15-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetCollagenMesh(MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>* pMesh_ab1)
{
    mpMesh_ab1 = pMesh_ab1;
	std::cout<<"hello2"<<std::endl;
}

//added 7-15-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetForceCollection_ab1(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces_ab1)
{
    //mpForceCollection = pForces;
	mpForceCollection_ab1 = pForces_ab1;
	//std::cout<<"pForces_ab->size()\t"<<pForces_ab->size()<<std::endl;
}


//added 9-13-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetCellForceOutputFile(out_stream pCellNodeForcesFile)
{
    mpCellNodeForcesFile = pCellNodeForcesFile;
}



//added 8-27-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetSubstrateForceOutputFile(out_stream pSubForcesFile)
{
    mpSubForcesFile = pSubForcesFile;
}



//added 6-22-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetSubstrateMesh(MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>* pMesh_ab)
{
    mpMesh_ab = pMesh_ab;
}

//added 8-4-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetForceCollection_ab(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces_ab)
{
    //mpForceCollection = pForces;
	mpForceCollection_ab = pForces_ab;
	//std::cout<<"pForces_ab->size()\t"<<pForces_ab->size()<<std::endl;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetForceCollection(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces)
{
    mpForceCollection = pForces;
	
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetUseAdaptiveTimestep(bool useAdaptiveTimestep)
{
    mUseAdaptiveTimestep = useAdaptiveTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::HasAdaptiveTimestep()
{
    return mUseAdaptiveTimestep;
}

//added 7-22-21

//*********************Computes forces for Collagen Mesh (control coming from UpdateCollagenTopology() in AbstractCellBasedSimulation/OffLatticeSimulation_ab1 and UpdateAllNodePositions_ab1() in ForwardEulerNumericalMethod_ab********************************************

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::ComputeForcesIncludingDamping_ab1()
{
    //CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

	//std::cout<<"inside ComputeForcesIncludingDamping_ab1()-cumulative forces for collagen in abstractnumericalmethod"<<std::endl;
	
	

	//exit(0);
	//std::cout<<"inside computeforcesincludingdamping_ab1"<<std::endl;		
	//exit(0);

    /*for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpMesh_ab1->GetNodeIteratorBegin();
         node_iter != mpMesh_ab1->GetNodeIteratorEnd(); ++node_iter)
    {
	   std::cout<<"inside clearappliedforce"<<std::endl;     
	node_iter->ClearAppliedForce();
    }*/
	//std::cout<<"loveandhate266220000000"<<std::endl;
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mpForceCollection_ab1->begin();
        iter != mpForceCollection_ab1->end(); ++iter)
    {
        //std::cout<<"inside mpforcecollection loop"<<std::endl;
        //std::cout<<iter<<std::endl;
	(*iter)->AddForceContribution_ab3(*mpMesh_ab1, *mpMesh_ab);   //added 6-25-20 
	   //std::cout<<"out"<<std::endl;
	//(*iter)->AddForceContribution(*mpCellPopulation);
    }
	//exit(0); 
    //std::cout<<"out"<<std::endl;
    /**
     * Here we deal with the special case forces on ghost nodes. Note that 'particles'
     * are dealt with like normal cells.
     *
     * \todo #2087 Consider removing dynamic_cast and instead adding a virtual method ApplyForcesToNonCellNodes()
     */
   /* if (mGhostNodeForcesEnabled)
    {
        dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(mpCellPopulation)->ApplyGhostForces();
    }*/
	//std::cout<<"out1"<<std::endl;
    // Store applied forces in a vector
    std::vector<c_vector<double, SPACE_DIM> > forces_as_vector;
    forces_as_vector.reserve(mpMesh_ab1->GetNumNodes());
	//std::cout<<"out2"<<std::endl;
	unsigned index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpMesh_ab1->GetNodeIteratorBegin();
         node_iter != mpMesh_ab1->GetNodeIteratorEnd(); ++node_iter, ++index)
    {
       

	double damping = mpMesh_ab1->GetDampingConstant();
	//double damping = 2;
        forces_as_vector.push_back(node_iter->rGetAppliedForce()/damping);
	//forces_as_vector.push_back(node_iter->rGetAppliedForce());

	

	//std::cout<<"totalforces_for_collagen[0]\t"<<forces_as_vector[index][0]<<"\t"<<"totalforces_for_collagen[1]\t"<<forces_as_vector[index][1]<<std::endl;

    }
	//std::cout<<"forces_as_vector.size()\t"<<forces_as_vector.size()<<std::endl;
   // CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

	
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpMesh_ab1->GetNodeIteratorBegin();
         node_iter != mpMesh_ab1->GetNodeIteratorEnd(); ++node_iter)
    	{
        	node_iter->ClearAppliedForce();
    	}



	//std::cout<<"inside computeforcesincludingdamping_ab1 - after force application"<<std::endl;
    return forces_as_vector;
	//std::cout<<"endofcomputeforcesincludingdamping@abstactnummeth"<<std::endl;
}

//*********************************************************************************************************************************************************


//added 8-3-20
//*********************Computes forces for Substrate Mesh (control coming from UpdateSubstrateTopology() in AbstractCellBasedSimulation/OffLatticeSimulation_ab1 and UpdateAllNodePositions_ab() in ForwardEulerNumericalMethod_ab********************************************

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::ComputeForcesIncludingDamping_ab()
{
    //CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);
	//std::cout<<"inside ComputeForcesIncludingDamping_ab1()-cumulative forces for substrate in abstractnumericalmethod"<<std::endl;
	//std::cout<<"ComputeForcesIncludingDamping_ab"<<std::endl;
	//exit(0);


	//std::string substrate_results_directory = "1-DVertexBasedMonolayer/substrate_results_from_time_";    //added 8-27-21

	//OutputFileHandler output_file_handler_ab3(substrate_results_directory+"/", true);    //added 8-27-21
	
	//mpSubForcesFile = output_file_handler_ab3.OpenOutputFile("substratenodeforces.dat");    //added 8-27-21

	//exit(0);

	//std::cout<<"b4_mpForceCollection_ab"<<std::endl;
	//std::cout<<"mpForceCollection_ab.size()\t"<<mpForceCollection_ab->size()<<std::endl;
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mpForceCollection_ab->begin();
        iter != mpForceCollection_ab->end(); ++iter)
    {
        //std::cout<<"inside mpforcecollection loop of ComputeForcesIncludingDamping_ab"<<std::endl;
        //std::cout<<iter<<std::endl;
	(*iter)->AddForceContribution_ab1(*mpMesh_ab);   //added 8-4-20 
	//(*iter)->AddForceContribution_ab(*mpCellPopulation, *mpMesh_ab);   //added 6-25-20     
	//(*iter)->AddForceContribution(*mpCellPopulation);
    }
    //std::cout<<"out"<<std::endl;
    /**
     * Here we deal with the special case forces on ghost nodes. Note that 'particles'
     * are dealt with like normal cells.
     *
     * \todo #2087 Consider removing dynamic_cast and instead adding a virtual method ApplyForcesToNonCellNodes()
     */
   /* if (mGhostNodeForcesEnabled)
    {
        dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(mpCellPopulation)->ApplyGhostForces();
    }*/

    // Store applied forces in a vector
    std::vector<c_vector<double, SPACE_DIM> > forces_as_vector_ab;
    forces_as_vector_ab.reserve(mpMesh_ab->GetNumNodes());


	unsigned index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpMesh_ab->GetNodeIteratorBegin();
         node_iter != mpMesh_ab->GetNodeIteratorEnd(); ++node_iter, ++index)
    {
        double damping = mpMesh_ab->GetDampingConstant();
	//double damping = 2;	
        forces_as_vector_ab.push_back(node_iter->rGetAppliedForce()/damping);
	//forces_as_vector_ab.push_back(node_iter->rGetAppliedForce());
	//std::cout<<"ComputeForcesIncludingDamping_ab2"<<std::endl;
	
	*mpSubForcesFile << forces_as_vector_ab[index][0] << " ";	//added 8-27-21		commented 3-25-22	

	std::cout<<"totalforces_for_substrate[0]\t"<<forces_as_vector_ab[index][0]<<"\t"<<"totalforces_for_substrate[1]\t"<<forces_as_vector_ab[index][1]<<std::endl;			//uncommented 10/19/22 for figuring out units
	//exit(0);

    }
	*mpSubForcesFile << "\n";	//added 8-27-21		commented 3-25-22
	//mpSubForcesFile->close(); 	//added 8-27-21
	//exit(0);

	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpMesh_ab->GetNodeIteratorBegin();
         node_iter != mpMesh_ab->GetNodeIteratorEnd(); ++node_iter)
    	{
        	node_iter->ClearAppliedForce();
    	}

    //CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

    return forces_as_vector_ab;
}

//*********************************************************************************************************************************************************



//*********************Computes forces for Cell Mesh (control coming from UpdateCellLocationAndTopology() in AbstractCellBasedSimulation********************************************
//modified 12-22-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::ComputeForcesIncludingDamping()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

	//std::cout<<"inside ComputeForcesIncludingDamping()-cumulative forces for cell in abstractnumericalmethod"<<std::endl;

    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != mpCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mpForceCollection->begin();
        iter != mpForceCollection->end(); ++iter)
    {
        //std::cout<<"loveandhate266220000000"<<std::endl;
        //std::cout<<iter<<std::endl;
	(*iter)->SetCellForceOutputFile(mpCellNodeForcesFile);   //added 10-6-21	

	(*iter)->AddForceContribution_ab2(*mpCellPopulation, *mpMesh_ab1);   //added 6-25-20 
	   
	//(*iter)->AddForceContribution(*mpCellPopulation);
    }
	//exit(0); 
    //std::cout<<"out"<<std::endl;
    /**
     * Here we deal with the special case forces on ghost nodes. Note that 'particles'
     * are dealt with like normal cells.
     *
     * \todo #2087 Consider removing dynamic_cast and instead adding a virtual method ApplyForcesToNonCellNodes()
     */
    if (mGhostNodeForcesEnabled)
    {
        dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(mpCellPopulation)->ApplyGhostForces();
    }
	//std::cout<<"out1"<<std::endl;
    // Store applied forces in a vector
    std::vector<c_vector<double, SPACE_DIM> > forces_as_vector;
    forces_as_vector.reserve(mpCellPopulation->GetNumNodes());
	//std::cout<<"out2"<<std::endl;

	unsigned index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != mpCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
        //double damping = mpCellPopulation->rGetMesh().GetDampingConstant();		//commented 12-22-21
	double damping = node_iter->rGetDampingConstant();		//added 12-22-21
	std::cout<<"damping\t"<<damping<<std::endl;	//added 10/10/22 for figuring out units	
	//double damping = 2;
        forces_as_vector.push_back(node_iter->rGetAppliedForce()/damping);

	//std::cout<<"damping for node is\t"<<damping<<std::endl;

	//*mpCellNodeForcesFile << forces_as_vector[index][0] << " ";		//added 9-13-21

	//std::cout<<"totalforces_for_cell[0]\t"<<forces_as_vector[index][0]<<"\t"<<"totalforces_for_cell[1]\t"<<forces_as_vector[index][1]<<std::endl;


    }
	//*mpCellNodeForcesFile << "\n";	//added 8-27-21

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);
	//std::cout<<"out3"<<std::endl;
    return forces_as_vector;
	//std::cout<<"endofcomputeforcesincludingdamping@abstactnummeth"<<std::endl;
}

//*********************************************************************************************************************************************************


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SaveCurrentLocations()
{
    std::vector<c_vector<double, SPACE_DIM> > current_locations;
    current_locations.reserve(mpCellPopulation->GetNumNodes());

    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        current_locations.push_back(node_iter->rGetLocation());
    }

    return current_locations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SafeNodePositionUpdate( unsigned nodeIndex, c_vector<double, SPACE_DIM> newPosition)
{
    ChastePoint<SPACE_DIM> new_point(newPosition);
    mpCellPopulation->SetNode(nodeIndex, new_point);
}

//added 8-8-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SafeNodePositionUpdate_ab( unsigned nodeIndex, c_vector<double, SPACE_DIM> newPosition)
{
    ChastePoint<SPACE_DIM> new_point(newPosition);
    mpMesh_ab->SetNode(nodeIndex, new_point);
}

//added 7-28-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SafeNodePositionUpdate_ab1( unsigned nodeIndex, c_vector<double, SPACE_DIM> newPosition)
{
    ChastePoint<SPACE_DIM> new_point(newPosition);
    mpMesh_ab1->SetNode(nodeIndex, new_point);
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::DetectStepSizeExceptions(unsigned nodeIndex, c_vector<double,SPACE_DIM>& displacement, double dt)
{
    try
    {
        mpCellPopulation->CheckForStepSizeException(nodeIndex, displacement, dt);
    }
    catch (StepSizeException& e)
    {
        if (!(e.IsTerminal()) && (mUseAdaptiveTimestep==false))
        {
            /*
             * If adaptivity is turned off but the simulation can continue, just produce a warning.
             * Only the case for vertex-based cell populations, which can alter node displacement directly
             * to avoid cell rearrangement problems.
             */
            WARN_ONCE_ONLY(e.what());
        }
        else
        {
            throw e;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetUseUpdateNodeLocation(bool useUpdateNodeLocation)
{
    mUseUpdateNodeLocation = useUpdateNodeLocation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::GetUseUpdateNodeLocation()
{
    return mUseUpdateNodeLocation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodInfo(out_stream& rParamsFile)
{
    std::string numerical_method_type = GetIdentifier();

    *rParamsFile << "\t\t<" << numerical_method_type << ">\n";
    OutputNumericalMethodParameters(rParamsFile);
    *rParamsFile << "\t\t</" << numerical_method_type << ">\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<UseAdaptiveTimestep>" << mUseAdaptiveTimestep << "</UseAdaptiveTimestep> \n";
    *rParamsFile << "\t\t\t<UseUpdateNodeLocation>" << mUseUpdateNodeLocation << "</UseUpdateNodeLocation> \n";
    *rParamsFile << "\t\t\t<GhostNodeForcesEnabled>" << mGhostNodeForcesEnabled << "</GhostNodeForcesEnabled> \n";
}

// Explicit instantiation
template class AbstractNumericalMethod<1,1>;
template class AbstractNumericalMethod<1,2>;
template class AbstractNumericalMethod<2,2>;
template class AbstractNumericalMethod<1,3>;
template class AbstractNumericalMethod<2,3>;
template class AbstractNumericalMethod<3,3>;

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

#include "ForwardEulerNumericalMethod_ab.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ForwardEulerNumericalMethod_ab<ELEMENT_DIM,SPACE_DIM>::ForwardEulerNumericalMethod_ab()
    : AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ForwardEulerNumericalMethod_ab<ELEMENT_DIM,SPACE_DIM>::~ForwardEulerNumericalMethod_ab()
{
}
/*
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetForceCollection(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces)
{
    mpForceCollection = pForces;
}

    void SetForceCollection(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces);

 std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* mpForceCollection;

*/
/*
//added 6-21-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethod_ab<ELEMENT_DIM,SPACE_DIM>::SetSubstrateMesh(MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>* pMesh_ab)
{
    mpMesh_ab = pMesh_ab;
}

*/


//added 7-22-21

//*********************Updates positions for Collagen Mesh (control coming from UpdateCollagenTopology() in AbstractCellBasedSimulation and OffLatticeSimulation_ab1. Control then goes to member function ComputeForcesIncludingDamping() of class AbstractNumericalMethod ********************************************

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethod_ab<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions_ab1(double dt)
{
    //std::cout<<"inside UpdateAllNodePositions_ab1 in forwardeulernumericalmethod_ab"<<std::endl;	
	//exit(0);
    //std::cout<<this<<std::endl;
    if (!this->mUseUpdateNodeLocation)
    {
        /*
	//{ added 6-22-20
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpMesh_ab->GetNodeIteratorBegin();
             node_iter != this->mpMesh_ab->GetNodeIteratorEnd();
             ++node_iter)
        {
		//std::cout<<"He-Manxyz"<<std::endl; 
		//std::vector<double>& abc = (node_iter)->rGetNodeAttributes();
		//std::cout<<abc[0]<<std::endl;           
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

	//}*/
        

        // Apply forces to each cell, and save a vector of net forces F
	//std::vector<c_vector<double, SPACE_DIM> > forces = this->ComputeForcesIncludingDamping();
        std::vector<c_vector<double, SPACE_DIM> > forces = this->ComputeForcesIncludingDamping_ab1();

        unsigned index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpMesh_ab1->GetNodeIteratorBegin();
             node_iter != this->mpMesh_ab1->GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            // Get the current node location and calculate the new location according to the forward Euler method
		//std::cout<<"UpdateAllNodePositions_ab1.5"<<std::endl;
            const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation();
		//std::cout<<"UpdateAllNodePositions_ab1.6"<<std::endl;
		//double abc = 10000000;
		c_vector<double, SPACE_DIM> displacement = dt * forces[index];
           // c_vector<double, SPACE_DIM> displacement = dt * forces[index]*abc;
		//std::cout<<"[index]\t"<<index<<"\tforces[0]\t"<<forces[index][0]<<std::endl;

            // In the vertex-based case, the displacement may be scaled if the cell rearrangement threshold is exceeded
            //this->DetectStepSizeExceptions(node_iter->GetIndex(), displacement, dt);
		//std::cout<<"UpdateAllNodePositions_ab1.8"<<std::endl;
            c_vector<double, SPACE_DIM> new_location = r_old_location + displacement;
            this->SafeNodePositionUpdate_ab1(node_iter->GetIndex(), new_location);
		//std::cout<<"UpdateAllNodePositions_ab2"<<std::endl;
        }
    }
	
}

//************************************************************************************************************************************************


//added 8-3-20

//*********************Updates positions for Substrate Mesh (control coming from UpdateSubstrateTopology() in AbstractCellBasedSimulation and OffLatticeSimulation_ab1. Control then goes to member function ComputeForcesIncludingDamping() of class AbstractNumericalMethod )********************************************

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethod_ab<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions_ab(double dt)
{
    //std::cout<<"UpdateAllNodePositions_ab"<<std::endl;	
    //std::cout<<this<<std::endl;
    if (!this->mUseUpdateNodeLocation)
    {
        /*
	//{ added 6-22-20
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpMesh_ab->GetNodeIteratorBegin();
             node_iter != this->mpMesh_ab->GetNodeIteratorEnd();
             ++node_iter)
        {
		//std::cout<<"He-Manxyz"<<std::endl; 
		//std::vector<double>& abc = (node_iter)->rGetNodeAttributes();
		//std::cout<<abc[0]<<std::endl;           
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

	//}*/
        

        // Apply forces to each cell, and save a vector of net forces F
	//std::vector<c_vector<double, SPACE_DIM> > forces = this->ComputeForcesIncludingDamping();
        std::vector<c_vector<double, SPACE_DIM> > forces = this->ComputeForcesIncludingDamping_ab();

        unsigned index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpMesh_ab->GetNodeIteratorBegin();
             node_iter != this->mpMesh_ab->GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            // Get the current node location and calculate the new location according to the forward Euler method
		//std::cout<<"UpdateAllNodePositions_ab1.5"<<std::endl;
            const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation();
		//std::cout<<"UpdateAllNodePositions_ab1.6"<<std::endl;
		//double abc = 10000000;
		c_vector<double, SPACE_DIM> displacement = dt * forces[index];
           // c_vector<double, SPACE_DIM> displacement = dt * forces[index]*abc;
	//	std::cout<<"displacement[0]\t"<<displacement[0]<<"\tdisplacement[1]\t"<<displacement[1]<<std::endl;

            // In the vertex-based case, the displacement may be scaled if the cell rearrangement threshold is exceeded
            //this->DetectStepSizeExceptions(node_iter->GetIndex(), displacement, dt);
		//std::cout<<"UpdateAllNodePositions_ab1.8"<<std::endl;
            c_vector<double, SPACE_DIM> new_location = r_old_location + displacement;
            this->SafeNodePositionUpdate_ab(node_iter->GetIndex(), new_location);
		//std::cout<<"UpdateAllNodePositions_ab2"<<std::endl;
        }
    }
	
}
//************************************************************************************************************************************************


//*********************Updates positions for Cell Mesh (control coming from UpdateCellLocationAndTopology() in AbstractCellBasedSimulation and OffLatticeSimulation_ab1. Control then goes to member function ComputeForcesIncludingDamping() of class AbstractNumericalMethod *******************************************

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethod_ab<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt)
{
    //std::cout<<"loveandhate2222"<<std::endl;	
    //std::cout<<this<<std::endl;
    if (!this->mUseUpdateNodeLocation)
    {
        /*
	//{ added 6-22-20
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpMesh_ab->GetNodeIteratorBegin();
             node_iter != this->mpMesh_ab->GetNodeIteratorEnd();
             ++node_iter)
        {
		//std::cout<<"He-Manxyz"<<std::endl; 
		//std::vector<double>& abc = (node_iter)->rGetNodeAttributes();
		//std::cout<<abc[0]<<std::endl;           
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

	//}*/
        

        // Apply forces to each cell, and save a vector of net forces F
        std::vector<c_vector<double, SPACE_DIM> > forces = this->ComputeForcesIncludingDamping();

        unsigned index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            // Get the current node location and calculate the new location according to the forward Euler method
            const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation();
            c_vector<double, SPACE_DIM> displacement = dt * forces[index];
		
		std::cout<<"dt\t"<<dt<<std::endl;	//added 10/10/22 for figuring out units
		std::cout<<"index\t"<<index<<std::endl;	//added 10/10/22 for figuring out units
		std::cout<<"forces[index]\t"<<forces[index][0]<<std::endl;	//added 10/10/22 for figuring out units
		std::cout<<"displacement\t"<<displacement[0]<<std::endl;	//added 10/10/22 for figuring out units
		//exit(0);	//added 10/10/22 for figuring out units
		
		double stop_disp = 0.0001;	//added 8-13-23, value before was 0.00001		

		//added 12-7-21

		if (node_iter->ProtrusiveForceOn()==true)
		{
			node_iter->SetUpdateAdhesion1(true);		//added 12-13-21 (once protrusion has happened, update adhesion1 will take place. After updateadhesion1 has happened, updateadhesion1 will be set to off)
			std::cout<<"displacement[0]\t"<<displacement[0]<<std::endl;
			//std::cout<<"node_iter->UpdateAdhesion1On()\t"<<node_iter->UpdateAdhesion1On()<<std::endl;	
		}

		if (node_iter->IsLeaderNode()==true && abs(displacement[0])<=stop_disp)	//modified on 8-13-23, intial value is 0.00001
		{
			//exit(0);
			node_iter->SetProtrusiveForce(true);
			//node_iter->SetUpdateAdhesion1(true);		//added 12-13-21
			std::cout<<"new_pos[0]\t"<<r_old_location[0]+displacement[0]<<std::endl;	
		}
		else if (node_iter->IsLeaderNode()==true && abs(displacement[0])>stop_disp)	////modified on 8-13-23, intial value is 0.00001
		{
			
			node_iter->SetProtrusiveForce(false);
			//node_iter->SetUpdateAdhesion1(false);		//added 12-13-21
			std::cout<<"new_pos[0]\t"<<r_old_location[0]+displacement[0]<<std::endl;
		}
		//std::cout<<"inside UpdateAllNodePositions()-in forwardeulernumericalmethod\t"<<"dt is\t"<<dt<<std::endl;
		
		//std::cout<<"displacement[0]\t"<<displacement[0]<<"\tdisplacement[1]\t"<<displacement[1]<<std::endl;

		//std::cout<<"forces[0]\t"<<forces[index][0]<<"\tforces[1]\t"<<forces[index][1]<<std::endl;

            // In the vertex-based case, the displacement may be scaled if the cell rearrangement threshold is exceeded
           // this->DetectStepSizeExceptions(node_iter->GetIndex(), displacement, dt);
		//std::cout<<"cell node index is@updateallnodepositiions@forwardeuler\t"<<(node_iter)->GetIndex()<<std::endl;
            c_vector<double, SPACE_DIM> new_location = r_old_location + displacement;
            this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
		//std::cout<<"UpdateAllNodePositions@mpcellpopulation"<<std::endl;
        }
    }
   /**/ 
}
//************************************************************************************************************************************************




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethod_ab<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(rParamsFile);
}

// Explicit instantiation
template class ForwardEulerNumericalMethod_ab<1,1>;
template class ForwardEulerNumericalMethod_ab<1,2>;
template class ForwardEulerNumericalMethod_ab<2,2>;
template class ForwardEulerNumericalMethod_ab<1,3>;
template class ForwardEulerNumericalMethod_ab<2,3>;
template class ForwardEulerNumericalMethod_ab<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ForwardEulerNumericalMethod_ab)

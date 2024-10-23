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

#include <cassert>

#include "Node.hpp"
#include "Exception.hpp"

//////////////////////////////////////////////////////////////////////////
// Constructors
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::CommonConstructor(unsigned index, bool isBoundaryNode)
{
    mIndex = index;
    mIsBoundaryNode = isBoundaryNode;
    mIsInternal = false;
    mIsDeleted = false;
    mpNodeAttributes = nullptr;
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, ChastePoint<SPACE_DIM> point, bool isBoundaryNode)
{
    mLocation = point.rGetLocation();
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, std::vector<double> coords, bool isBoundaryNode)
{
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        mLocation(i) = coords.at(i);
    }
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, c_vector<double, SPACE_DIM> location, bool isBoundaryNode)
{
    mLocation = location;
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, bool isBoundaryNode, double v1, double v2, double v3)
{
    mLocation[0] = v1;
    if (SPACE_DIM > 1)
    {
        mLocation[1] = v2;
        if (SPACE_DIM > 2)
        {
            mLocation[2] = v3;
        }
    }
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, double *location, bool isBoundaryNode)
{
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        mLocation(i) = location[i];
    }
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::~Node()
{
    delete mpNodeAttributes;
}

//////////////////////////////////////////////////////////////////////////
// Methods dealing with node location
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetPoint(ChastePoint<SPACE_DIM> point)
{
    mLocation = point.rGetLocation();
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetIndex(unsigned index)
{
    mIndex = index;
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetAsBoundaryNode(bool value)
{
    mIsBoundaryNode = value;
}

template<unsigned SPACE_DIM>
ChastePoint<SPACE_DIM> Node<SPACE_DIM>::GetPoint() const
{
    return ChastePoint<SPACE_DIM>(mLocation);
}

template<unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& Node<SPACE_DIM>::rGetLocation() const
{
    // This assert statement is a useful warning: when new nodes are created we overwrite previously deleted nodes if there are any.
    // This means that we can not use this method to interrogate deleted nodes about their position before deletion because we can't
    // guarantee that the node has not been overwritten already. Hence, when implementing new functionality we need to make sure
    // that this functionality does not rely on being able to interrogate deleted nodes for their location.
    // \todo #2401: make this an exception.
    assert(!mIsDeleted);
    return mLocation;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM>& Node<SPACE_DIM>::rGetModifiableLocation()
{
    assert(!mIsDeleted);
    return mLocation;
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetIndex() const
{
    return mIndex;
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsBoundaryNode() const
{
    return mIsBoundaryNode;
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddNodeAttribute(double attribute)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddAttribute(attribute);
}

template<unsigned SPACE_DIM>
std::vector<double>& Node<SPACE_DIM>::rGetNodeAttributes()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetAttributes();
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetNumNodeAttributes()
{
    unsigned num_attributes;
    if (!mpNodeAttributes)
    {
        num_attributes = 0u;
    }
    else
    {
        num_attributes = mpNodeAttributes->rGetAttributes().size();
    }

    return num_attributes;
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::HasNodeAttributes()
{
    return (mpNodeAttributes != nullptr);
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM>& Node<SPACE_DIM>::rGetAppliedForce()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetAppliedForce();
}


//added 7-27-21   checking existance of bond between collagen and substrate
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::HasFocal1(bool flag2)
{
    ConstructNodeAttributes();

    mpNodeAttributes->HasFocal1(flag2);
}


//added 12-13-21   setting attribute of the node to be able to updateadhesion1 or otherwise 
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetUpdateAdhesion1(bool flag5)
{
    ConstructNodeAttributes();

    mpNodeAttributes->SetUpdateAdhesion1(flag5);
}



//added 12-7-21   setting attribute of the node to having protrusive force or otherwise 
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetProtrusiveForce(bool flag4)
{
    ConstructNodeAttributes();

    mpNodeAttributes->SetProtrusiveForce(flag4);
}


//added 9-16-21   setting attribute of the node to leader cell node or otherwise 
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::HasLeaderNode(bool flag3)
{
    ConstructNodeAttributes();

    mpNodeAttributes->HasLeaderNode(flag3);
}

//added 7-27-21   checking existance of bond between collagen and substrate
template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsFocal1()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->IsFocal1();
}


//added 12-13-21   checking whether UpdateAdhesion1 should be applied or not
template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::UpdateAdhesion1On()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->UpdateAdhesion1On();
}



//added 12-7-21   checking whether protrusive force should be applied or not
template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::ProtrusiveForceOn()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->ProtrusiveForceOn();
}



//added 9-16-21   checking whether node is a leader cell node or not
template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsLeaderNode()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->IsLeaderNode();
}





//added 7-14-20
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::HasFocal(bool flag1)
{
    ConstructNodeAttributes();

    mpNodeAttributes->HasFocal(flag1);
}

//added 7-14-20
template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsFocal()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->IsFocal();
}



//added 7-27-21
template<unsigned SPACE_DIM>
unsigned& Node<SPACE_DIM>::rGetCellNodeIndex()
{
	//std::cout<<"inside_rgetcollagennodeindex_beforeattributecheck"<<std::endl;    
	CheckForNodeAttributes();
	//std::cout<<"inside_rgetcollagennodeindex_afterattributecheck"<<std::endl;
	//std::cout<<"mpNodeAttributes->rGetCollagenNodeIndex()\t"<<mpNodeAttributes->rGetCollagenNodeIndex()<<std::endl;
    return mpNodeAttributes->rGetCellNodeIndex();
	
	
}


//added 7-27-21
template<unsigned SPACE_DIM>
double& Node<SPACE_DIM>::rGetCellNodeDistance()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetCellNodeDistance();
}


//added 7-27-21
template<unsigned SPACE_DIM>
std::vector<double>& Node<SPACE_DIM>::rGetCellNodeLocation()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetCellNodeLocation();
}






//added 7-16-21
template<unsigned SPACE_DIM>
unsigned& Node<SPACE_DIM>::rGetCollagenNodeIndex()
{
	//std::cout<<"inside_rgetcollagennodeindex_beforeattributecheck"<<std::endl;    
	CheckForNodeAttributes();
	//std::cout<<"inside_rgetcollagennodeindex_afterattributecheck"<<std::endl;
	//std::cout<<"mpNodeAttributes->rGetCollagenNodeIndex()\t"<<mpNodeAttributes->rGetCollagenNodeIndex()<<std::endl;
    return mpNodeAttributes->rGetCollagenNodeIndex();
	
	
}


//added 12-22-21
template<unsigned SPACE_DIM>
double& Node<SPACE_DIM>::rGetDampingConstant()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetDampingConstant();
}



//added 11-26-21
template<unsigned SPACE_DIM>
double& Node<SPACE_DIM>::rGetElementStiffnessContribution()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetElementStiffnessContribution();
}



//added 9-20-21
template<unsigned SPACE_DIM>
double& Node<SPACE_DIM>::rGetCollagenCellRestLength()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetCollagenCellRestLength();
}


//added 7-16-21
template<unsigned SPACE_DIM>
double& Node<SPACE_DIM>::rGetCollagenNodeDistance()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetCollagenNodeDistance();
}


//added 7-16-21
template<unsigned SPACE_DIM>
std::vector<double>& Node<SPACE_DIM>::rGetCollagenNodeLocation()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetCollagenNodeLocation();
}






//added 7-16-20
template<unsigned SPACE_DIM>
unsigned& Node<SPACE_DIM>::rGetSubstrateNodeIndex()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetSubstrateNodeIndex();
}




//added 7-14-20
template<unsigned SPACE_DIM>
double& Node<SPACE_DIM>::rGetSubstrateNodeDistance()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetSubstrateNodeDistance();
}



//added 7-5-20
template<unsigned SPACE_DIM>
std::vector<double>& Node<SPACE_DIM>::rGetSubstrateNodeLocation()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetSubstrateNodeLocation();
}

//added 7-7-20
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::ClearNodeLocation()
{
    ConstructNodeAttributes();

    mpNodeAttributes->ClearNodeLocation();
}


template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::ClearAppliedForce()
{
    ConstructNodeAttributes();

    mpNodeAttributes->ClearAppliedForce();
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddAppliedForceContribution(const c_vector<double, SPACE_DIM>& rForceContribution)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddAppliedForceContribution(rForceContribution);
}


//added 12-22-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddDampingConstant(double& rDampingConstant)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddDampingConstant(rDampingConstant);
}



//added 11-26-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddElementStiffnessContribution(double& rElementStiffnessContribution)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddElementStiffnessContribution(rElementStiffnessContribution);
}


//added 7-27-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddCellNodeLocation(const c_vector<double, SPACE_DIM>& rCellNodeLocation)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddCellNodeLocation(rCellNodeLocation);
}

//added 7-27-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddCellNodeDistance(double& rCellNodeDistance)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddCellNodeDistance(rCellNodeDistance);
}

//added 7-27-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddCellNodeIndex(unsigned& rCellNodeIndex)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddCellNodeIndex(rCellNodeIndex);
}




//added 7-15-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddCollagenNodeLocation(const c_vector<double, SPACE_DIM>& rCollagenNodeLocation)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddCollagenNodeLocation(rCollagenNodeLocation);
}

//added 9-20-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddCollagenCellRestLength(double& rCollagenCellRestLength)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddCollagenCellRestLength(rCollagenCellRestLength);
}



//added 7-15-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddCollagenNodeDistance(double& rCollagenNodeDistance)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddCollagenNodeDistance(rCollagenNodeDistance);
}

//added 7-15-21
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddCollagenNodeIndex(unsigned& rCollagenNodeIndex)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddCollagenNodeIndex(rCollagenNodeIndex);
}





//added 7-5-20
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddSubstrateNodeLocation(const c_vector<double, SPACE_DIM>& rSubstrateNodeLocation)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddSubstrateNodeLocation(rSubstrateNodeLocation);
}

//added 7-14-20
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddSubstrateNodeDistance(double& rSubstrateNodeDistance)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddSubstrateNodeDistance(rSubstrateNodeDistance);
}


//added 7-14-20
template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddSubstrateNodeIndex(unsigned& rSubstrateNodeIndex)
{
    ConstructNodeAttributes();

    mpNodeAttributes->AddSubstrateNodeIndex(rSubstrateNodeIndex);
}



template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsParticle()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->IsParticle();
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetIsParticle(bool isParticle)
{
    ConstructNodeAttributes();

    mpNodeAttributes->SetIsParticle(isParticle);
}

template<unsigned SPACE_DIM>
double Node<SPACE_DIM>::GetRadius()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->GetRadius();
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetRadius(double radius)
{
    ConstructNodeAttributes();

    mpNodeAttributes->SetRadius(radius);
}

//////////////////////////////////////////////////////////////////////////
// Tracking (boundary) elements which contain this node as a vertex
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddElement(unsigned index)
{
    mElementIndices.insert(index);
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::RemoveElement(unsigned index)
{
    unsigned count = mElementIndices.erase(index);
    if (count == 0)
    {
        EXCEPTION("Tried to remove an index which was not in the set");
    }
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::RemoveBoundaryElement(unsigned index)
{
    unsigned count = mBoundaryElementIndices.erase(index);
    if (count == 0)
    {
        EXCEPTION("Tried to remove an index which was not in the set");
    }
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddBoundaryElement(unsigned index)
{
    mBoundaryElementIndices.insert(index);
}

template<unsigned SPACE_DIM>
std::set<unsigned>& Node<SPACE_DIM>::rGetContainingElementIndices()
{
    return mElementIndices;
}

template<unsigned SPACE_DIM>
std::set<unsigned>& Node<SPACE_DIM>::rGetContainingBoundaryElementIndices()
{
    return mBoundaryElementIndices;
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetNumContainingElements() const
{
    return mElementIndices.size();
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetNumBoundaryElements() const
{
    return mBoundaryElementIndices.size();
}

//////////////////////////////////////////////////////////////////////////
// Tracking neighbours of the node
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddNeighbour(unsigned index)
{
    ConstructNodeAttributes();

    return mpNodeAttributes->AddNeighbour(index);
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::ClearNeighbours()
{
    ConstructNodeAttributes();

    mpNodeAttributes->ClearNeighbours();
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::RemoveDuplicateNeighbours()
{
    CheckForNodeAttributes();

    mpNodeAttributes->RemoveDuplicateNeighbours();
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::NeighboursIsEmpty()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->NeighboursIsEmpty();
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetNeighboursSetUp(bool flag)
{
    ConstructNodeAttributes();

    mpNodeAttributes->SetNeighboursSetUp(flag);
};

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::GetNeighboursSetUp()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->GetNeighboursSetUp();
};

template<unsigned SPACE_DIM>
std::vector<unsigned>& Node<SPACE_DIM>::rGetNeighbours()
{
    CheckForNodeAttributes();

    return mpNodeAttributes->rGetNeighbours();
};

//////////////////////////////////////////////////////////////////////////
// Methods dealing with some node flags (deleted, region)
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::CheckForNodeAttributes() const
{
    if (mpNodeAttributes == nullptr)
    {
        EXCEPTION("Node has no attributes associated with it. Construct attributes first");
    }
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::ConstructNodeAttributes()
{
    if (mpNodeAttributes == nullptr)
    {
        mpNodeAttributes = new NodeAttributes<SPACE_DIM>();
    }
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::MarkAsDeleted()
{
    mIsDeleted = true;
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsDeleted() const
{
    return mIsDeleted;
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::MarkAsInternal()
{
    mIsInternal = true;
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsInternal() const
{
    return mIsInternal;
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetRegion(unsigned region)
{
    ConstructNodeAttributes();
    mpNodeAttributes->SetRegion(region);
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetRegion() const
{
    unsigned region = 0;

    if (mpNodeAttributes)
    {
        region = mpNodeAttributes->GetRegion();
    }

    return region;
}

// Explicit instantiation
template class Node<1>;
template class Node<2>;
template class Node<3>;

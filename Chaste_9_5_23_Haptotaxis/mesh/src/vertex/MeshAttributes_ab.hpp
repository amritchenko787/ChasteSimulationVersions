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

#ifndef _MESHATTRIBUTES_AB_HPP_
#define _MESHATTRIBUTES_AB_HPP_

#include "UblasVectorInclude.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
//#include "MutableVertexMesh.hpp"
/**
 * A container for attributes associated with the Node class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MeshAttributes_ab
{
private:

    /** Arbitrary attributes that a user gives meaning to */
    std::vector<double> mAttributes;

    /** The ID of the region of mesh in which the Node lies */
    unsigned mRegion;

    /** For mutable nodes in OffLatticeSimulations, a container for the force accumulated on this node. */
    c_vector<double, SPACE_DIM> mAppliedForce;
    
	std::vector<double> mSubstrateNodeLocation;  //added 7-5-20
	double mSubstrateNodeDistance;  //added 7-14-20
	unsigned mSubstrateNodeIndex;  //added 7-16-20


    /** The radius associated with the Node */
    double mRadius;

    /** Vector of indices corresponding to neighbouring nodes. */
    std::vector<unsigned> mNeighbourIndices;

    /** A bool indicating whether the neighbours of this node have been calculated yet. */
    bool mNeighboursSetUp;

    /** Whether the node represents a particle or not: Used for NodeBasedCellPopulationWithParticles */
    bool mIsParticle;

	bool mflag1;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mAttributes;
        archive & mRegion;
        archive & mRadius;
        archive & mNeighbourIndices;
        archive & mNeighboursSetUp;
        archive & mIsParticle;

        for (unsigned d = 0; d < SPACE_DIM; d++)
        {
            archive & mAppliedForce[d];
		archive & mSubstrateNodeLocation[d];     //added 7-5-20
        } 
    }

public:

    /**
     * Defaults all variables.
     */
    MeshAttributes_ab();

    /**
     * @return mAttributes
     */
    std::vector<double>& rGetAttributes();

    /**
     * Push an attribute back onto mAttributes
     *
     * @param attribute the value of the attribute.
     */
    void AddAttribute(double attribute);

    /**
     * Get the region ID
     *
     * @return mRegion
     */
    unsigned GetRegion();

    /**
     * Set the region ID
     *
     * @param region the value to to assign to mRegion.
     */
    void SetRegion(unsigned region);

    /**
     * Get the current value of the applied force on the node.
     *
     * @return mAppliedForce
     */
    c_vector<double, SPACE_DIM>& rGetAppliedForce();
	std::vector<double>& rGetSubstrateNodeLocation();    //added 7-5-20
	double& rGetSubstrateNodeDistance();    //added 7-14-20
	bool IsFocal();     //added 7-14-20
	unsigned& rGetSubstrateNodeIndex();    //added 7-16-20


    /**
     * Add a contribution to the force vector
     *
     * @param rForceContribution the contribution to add to mAppliedForce
     */
    void AddAppliedForceContribution(const c_vector<double, SPACE_DIM>& rForceContribution);
    void AddSubstrateNodeLocation(const c_vector<double, SPACE_DIM>& rSubstrateNodeLocation); //added 7-5-20
    void AddSubstrateNodeDistance(double& rSubstrateNodeDistance); //added 7-14-20
	void HasFocal(bool flag1);    //added 7-14-20
	void AddSubstrateNodeIndex(unsigned& rSubstrateNodeIndex); //added 7-16-20



    /**
     * Set mAppliedForce to a zero vector.
     */
    void ClearAppliedForce();
    void ClearNodeLocation();  //added 7-7-20
    /**
     * Add a neighbour to this node's vector of neighbouring node  indices.
     *
     * @param index of the node to add.
     */
    void AddNeighbour(unsigned index);

    /**
     * Clear this node's vector of neighbour indices.
     */
    void ClearNeighbours();

    /**
     * Remove duplicates from the vector of node neighbour indices.
     */
    void RemoveDuplicateNeighbours();

    /**
     * Check whether the node neighbours collection is empty.
     *
     * @return whether this node has any neighbours.
     */
    bool NeighboursIsEmpty();

    /**
     * Sets a flag to indicate that the neighbours of this node have/have not been updated.
     *
     * @param flag whether the neighbours are set up or not
     */
    void SetNeighboursSetUp(bool flag);

    /**
     * @return a flag to indicate that the neighbours of this node have/have not been updated.
     */
    bool GetNeighboursSetUp();

    /**
     * @return this node's vector of neighbour indices.
     */
    std::vector<unsigned>& rGetNeighbours();

    /**
     * Get whether this node is a particle, or not.
     *
     * @return mIsParticle
     */
    bool IsParticle();

    /**
     * Set the flag mIsParticle.
     * @param isParticle whether this node is particle or not.
     */
    void SetIsParticle(bool isParticle);

    /**
     * Return the radius associated with the Node
     *
     * @return mRadius
     */
    double GetRadius();

    /**
     * Set the value of the radius.
     *
     * @param radius the value to assign to mRadius. Must be >= 0.0
     */
    void SetRadius(double radius);
};

#endif //_NODEATTRIBUTES_HPP_
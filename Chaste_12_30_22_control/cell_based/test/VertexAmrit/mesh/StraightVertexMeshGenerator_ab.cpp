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

#include "StraightVertexMeshGenerator_ab.hpp"

StraightVertexMeshGenerator_ab::StraightVertexMeshGenerator_ab(unsigned numElementsAcross,
                                                           unsigned numElementsUp,
							   double ycoordinate,     //added 7-12-21	
							   bool isCollagen,     //added 7-12-21	
                                                           bool isFlatBottom,
                                                           double cellRearrangementThreshold,
                                                           double t2Threshold,
                                                           double elementArea)
{
    assert(numElementsAcross > 0);
	assert(numElementsAcross % 2 == 0);    //added 7-1-21
    //assert(numElementsUp > 0);    //mod 7-1-21
    //assert(cellRearrangementThreshold > 0.0);    //mod 7-1-21
    //assert(t2Threshold > 0.0);    //mod 7-1-21
    //assert(elementArea > 0.0);    //mod 7-1-21
    //std::cout<<"damnation"<<std::endl;
    
    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;
	std::vector<Node<2>*> element_nodes1;    //added 8-9-20
	std::vector<Node<2>*> element_nodes2;    //added 8-9-20
	//std::vector<Node<2>*> element_nodes3;    //added 8-9-20    //mod 7-1-21
	//std::vector<Node<2>*> element_nodes4;    //added 8-9-20    //mod 7-1-21
	//std::vector<Node<2>*> element_nodes5;    //added 8-9-20    //mod 7-1-21
	//std::vector<Node<2>*> element_nodes6;    //added 8-9-20    //mod 7-1-21

    unsigned node_index = 0;
    //unsigned node_indices[18];    //mod 7-1-21   
	unsigned node_indices[2];    //added 7-1-21
    unsigned element_index;
    //std::cout<<"damnation_2"<<std::endl;
    // Create the nodes, row by row, from the bottom up
	
	double kk = 0;
	double numelementsacross = numElementsAcross;
    // On the first row we have numElementsAcross nodes, all of which are boundary nodes
    for (unsigned i=0; i<numElementsAcross + 1; i++)
    {
       
	Node<2>* p_node = new Node<2>(node_index, true, kk - (numelementsacross/4), ycoordinate);
	//std::cout<<"damnation_abcd"<<"\t"<<kk<<"\t"<<kk - (numelementsacross/4)<<std::endl;
	//double abcd = 2;
	//p_node->SetRadius(abcd);
	//std::cout<<"damnation_abcd2"<<std::endl;
	p_node->HasFocal(false);
	
	p_node->HasFocal1(false);	//added 7-27-21

	//p_node->HasFocal(false);
        nodes.push_back(p_node);
        node_index++;
	kk = kk + 0.5;
    }

//std::cout<<"damnation_3"<<std::endl;

    /*
     * On each interior row we have numElementsAcross+1 nodes. On the second and penultimate
     * row all nodes are boundary nodes. On other rows the first and last nodes only
     * are boundary nodes.
     */

  

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    
	for (unsigned i=0; i<numElementsAcross; i++)
        	{
			node_indices[0] = i;
			node_indices[1] = i + 1;
			for (unsigned k=0; k<2; k++)
			{
				element_nodes1.push_back(nodes[node_indices[k]]);
			}
			element_index = i;
			elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
			element_nodes1.clear();			
			
		}
			



    mpMesh = new MutableVertexMesh<2,2>(nodes, elements, cellRearrangementThreshold, t2Threshold);

    // Scale the mesh so that each element's area takes the value elementArea
    mpMesh->Scale(sqrt(elementArea*2.0/sqrt(3.0)), sqrt(elementArea*2.0/sqrt(3.0)));
}

StraightVertexMeshGenerator_ab::~StraightVertexMeshGenerator_ab()
{
    delete mpMesh;
}

MutableVertexMesh<2,2>* StraightVertexMeshGenerator_ab::GetMesh()
{
    return mpMesh;
}

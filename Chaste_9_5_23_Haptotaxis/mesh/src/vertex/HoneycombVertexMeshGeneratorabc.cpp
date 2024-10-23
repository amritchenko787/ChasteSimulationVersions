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

#include "HoneycombVertexMeshGeneratorabc.hpp"

HoneycombVertexMeshGeneratorabc::HoneycombVertexMeshGeneratorabc(unsigned numElementsAcross,
                                                           unsigned numElementsUp,
                                                           bool isFlatBottom,
                                                           double cellRearrangementThreshold,
                                                           double t2Threshold,
                                                           double elementArea)
{
    assert(numElementsAcross > 0);
    assert(numElementsUp > 0);
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);
    assert(elementArea > 0.0);
    std::cout<<"damnation"<<std::endl;
    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;

    unsigned node_index = 0;
    unsigned node_indices[18];
    unsigned element_index;

    // Create the nodes, row by row, from the bottom up

    // On the first row we have numElementsAcross nodes, all of which are boundary nodes
    for (unsigned i=0; i<numElementsAcross; i++)
    {
       
	Node<2>* p_node = new Node<2>(node_index, true, i+0.5, 0);
        nodes.push_back(p_node);
        node_index++;
    }

    /*
     * On each interior row we have numElementsAcross+1 nodes. On the second and penultimate
     * row all nodes are boundary nodes. On other rows the first and last nodes only
     * are boundary nodes.
     */
    
    for (unsigned j=1; j<3*numElementsUp+1; j++)
    {
	if ((j%6 == 2)||(j%6 == 5))
        {
	        for (unsigned i=0; i<numElementsAcross; i++)
	    	{
                	double x_coord = (j%6 == 2) ? i+0.5 : i+1;
			double y_coord = j*0.5/sqrt(3.0);			    	  	 	 
			Node<2>* p_node = new Node<2>(node_index, false, x_coord, y_coord);
            		nodes.push_back(p_node);	
             		node_index++;
	    	}
	}
/*
	else if (j%6 == 5)
	{
		for (unsigned i=0; i<numElementsAcross: i++)
	    	{
                	double x_coord = i+1;
			double y_coord = j*0.5/sqrt(3.0);			    	  	 	 
			Node<2>* p_node = new Node<2>(node_index, false, x_coord, y_coord);
            		nodes.push_back(p_node);	
             		node_index++;
	    	}
	}	
*/	 
	else
	{      
		for (unsigned i=0; i<=numElementsAcross; i++)
        	{
            		double x_coord = ((j%6 == 0)||(j%6 == 4)) ? i+0.5 : i;
            		double y_coord = j*0.5/sqrt(3.0);
            		bool is_boundary_node = (j==1 || j==3*numElementsUp || i==0 || i==numElementsAcross) ? true : false;

            		Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord, y_coord);
            		nodes.push_back(p_node);
            		node_index++;
        	}
    	}
    }

/*
    for (unsigned j=1; j<2*numElementsUp+1; j++)
    {
        for (unsigned i=0; i<=numElementsAcross; i++)
        {
            double x_coord = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
            double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3.0);
            bool is_boundary_node = (j==1 || j==2*numElementsUp || i==0 || i==numElementsAcross) ? true : false;

            Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord, y_coord);
            nodes.push_back(p_node);
            node_index++;
        }
    }
*/
    /*
     * On the last row we have numElementsAcross nodes, all of which are boundary nodes.
     */
    double y_coord = (1.5*(2*numElementsUp+1) - 0.5*((2*numElementsUp+1)%2))*0.5/sqrt(3.0);
    if (((2*numElementsUp+1)%4 == 0)||((2*numElementsUp+1)%4 == 3))  //if NU is odd
    {
        Node<2>* p_node = new Node<2>(node_index, true, 0.5, y_coord);
        nodes.push_back(p_node);
        node_index++;
    }
    for (unsigned i=1; i<numElementsAcross; i++)
    {
        double x_coord = (((2*numElementsUp+1)%4 == 0)||((2*numElementsUp+1)%4 == 3)) ? i+0.5 : i;

        Node<2>* p_node = new Node<2>(node_index, true, x_coord, y_coord);
        nodes.push_back(p_node);
        node_index++;
    }
    if (((2*numElementsUp+1)%4 == 1)||((2*numElementsUp+1)%4 == 2))  //if NU is even
    {
        Node<2>* p_node = new Node<2>(node_index, true, numElementsAcross, y_coord);
        nodes.push_back(p_node);
        node_index++;
    }

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    
	for (unsigned j=0; j<numElementsUp; j++)
    	{
        	for (unsigned i=0; i<numElementsAcross; i++)
        	{
                	node_indices[0] = j*(5 + 3*(numElementsAcross - 1)) + i - 1*(j%2==0);
            		node_indices[1] = node_indices[0] + numElementsAcross + 1 + 1*(j%2==0 && j>0);
			node_indices[2] = node_indices[1] + numElementsAcross;
			
			std::vector<Node<2>*> element_nodes1;	
		        for (unsigned k=0; k<3; k++)
            		{
               			element_nodes1.push_back(nodes[node_indices[k]]);
            		}
			
			element_index = 6*j*numElementsAcross + 6*i + 1;			
			VertexElement<2,2>* p_element1 = new VertexElement<2,2>(element_index, element_nodes1);
            		elements.push_back(p_element1);
			element_index++;

			node_indices[3] = node_indices[1];
			node_indices[4] = node_indices[3] + 2*(numElementsAcross) + 1;
			node_indices[5] = node_indices[2];

			std::vector<Node<2>*> element_nodes2;	
		        for (unsigned k=3; k<6; k++)
            		{
               			element_nodes2.push_back(nodes[node_indices[k]]);
            		}
			VertexElement<2,2>* p_element2 = new VertexElement<2,2>(element_index, element_nodes2);
            		elements.push_back(p_element2);
			element_index++;

			node_indices[6] = node_indices[4];
			node_indices[7] = node_indices[6] + numElementsAcross + 1*(j%2 == 1);
			node_indices[8] = node_indices[2];

			std::vector<Node<2>*> element_nodes3;	
		        for (unsigned k=6; k<9; k++)
            		{
               			element_nodes3.push_back(nodes[node_indices[k]]);
            		}
			VertexElement<2,2>* p_element3 = new VertexElement<2,2>(element_index, element_nodes3);
            		elements.push_back(p_element3);
			element_index++;

			node_indices[9] = node_indices[7];
			node_indices[10] = node_indices[9] - (numElementsAcross + 1) - 1*(j%2 == 1);
			node_indices[11] = node_indices[2];
			
			std::vector<Node<2>*> element_nodes4;	
		        for (unsigned k=9; k<12; k++)
            		{
               			element_nodes4.push_back(nodes[node_indices[k]]);
            		}
			VertexElement<2,2>* p_element4 = new VertexElement<2,2>(element_index, element_nodes4);
            		elements.push_back(p_element4);
			element_index++;

			node_indices[12] = node_indices[10];
			node_indices[13] = node_indices[12] - 2*(numElementsAcross) - 1;
			node_indices[14] = node_indices[2];

			std::vector<Node<2>*> element_nodes5;	
		        for (unsigned k=12; k<15; k++)
            		{
               			element_nodes5.push_back(nodes[node_indices[k]]);
            		}
			VertexElement<2,2>* p_element5 = new VertexElement<2,2>(element_index, element_nodes5);
            		elements.push_back(p_element5);
			element_index++;

			node_indices[15] = node_indices[13];	
			node_indices[16] = node_indices[0];
			node_indices[17] = node_indices[2];

			std::vector<Node<2>*> element_nodes6;	
		        for (unsigned k=15; k<18; k++)
            		{
               			element_nodes6.push_back(nodes[node_indices[k]]);
            		}
			VertexElement<2,2>* p_element6 = new VertexElement<2,2>(element_index, element_nodes6);
            		elements.push_back(p_element6);
		}
	}




/*
    for (unsigned j=0; j<numElementsUp; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            if (j==0)
            {
                node_indices[0] = i;
            }
            else
            {
                node_indices[0] = 2*j*(numElementsAcross+1) - 1*(j%2==0) + i; // different for even/odd rows
            }
            node_indices[1] = node_indices[0] + numElementsAcross + 1 + 1*(j%2==0 && j>0);
            node_indices[2] = node_indices[1] + numElementsAcross + 1;
            node_indices[3] = node_indices[2] + numElementsAcross + 1*(j%2==1 && j<numElementsUp-1);
            node_indices[4] = node_indices[2] - 1;
            node_indices[5] = node_indices[1] - 1;

            std::vector<Node<2>*> element_nodes;
            for (unsigned k=0; k<6; k++)
            {
               element_nodes.push_back(nodes[node_indices[k]]);
            }

            element_index = j*numElementsAcross + i;
            VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
            elements.push_back(p_element);
        }
    }
*/

    mpMesh = new MutableVertexMesh<2,2>(nodes, elements, cellRearrangementThreshold, t2Threshold);

    // Scale the mesh so that each element's area takes the value elementArea
    mpMesh->Scale(sqrt(elementArea*2.0/sqrt(3.0)), sqrt(elementArea*2.0/sqrt(3.0)));
}

HoneycombVertexMeshGeneratorabc::~HoneycombVertexMeshGeneratorabc()
{
    delete mpMesh;
}

MutableVertexMesh<2,2>* HoneycombVertexMeshGeneratorabc::GetMesh()
{
    return mpMesh;
}

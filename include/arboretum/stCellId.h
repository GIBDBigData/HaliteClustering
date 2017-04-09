/**********************************************************************
 * GBDI Arboretum - Copyright (c) 2002-2009 GBDI-ICMC-USP
 *
 *                           Homepage: http://gbdi.icmc.usp.br/arboretum
 **********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2009 Grupo de Bases de Dados e Imagens, Instituto de
 * Ci�ncias Matem�ticas e de Computa���o, University of S�o Paulo -
 * Brazil (the Databases and Image Group - Intitute of Matematical and
 * Computer Sciences).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *    if any, must include the following acknowledgment:
 *       "This product includes software developed by Grupo de Bases
 *        de Dados e Imagens, Instituto de Ci�ncias Matem�ticas e de
 *        Computa��o, University of S�o Paulo - Brazil (the Databases 
 *        and Image Group - Intitute of Matematical and Computer
 *        Sciences)"
 *
 *    Alternately, this acknowledgment may appear in the software itself,
 *    if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names of the research group, institute, university, authors
 *    and collaborators must not be used to endorse or promote products
 *    derived from this software without prior written permission.
 *
 * 5. The names of products derived from this software may not contain
 *    the name of research group, institute or university, without prior
 *    written permission of the authors of this software.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OF THIS SOFTWARE OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * ====================================================================
 *                                            http://gbdi.icmc.usp.br/
 */
/**
 * @file
 * This file defines the class stCellId.
 *
 * @version 3.0
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
 * @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 * 
 */
// Copyright (c) 2002-2009 GBDI-ICMC-USP

#ifndef __STCELLID_H
#define __STCELLID_H

#include "stCommon.h"
#include <vector>
#include <memory>
#include <algorithm>

//----------------------------------------------------------------------------
// class stCellId
//----------------------------------------------------------------------------
/**
 * To identify cells in datasets.
 *
 * @version 3.0
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
 * @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 */
class stCellId {   

    public:

        /**
         * Vector that is going to be the index in which each bit refers to a dimension.
         * Possible unused bits will be in the "left side" of index[0].
         */
        std::vector<unsigned char> index; 

    public:

        /**
         * Constructor
         */
	stCellId() {
	}
        stCellId(const size_t DIM){         
            const size_t size = (DIM + 7) / 8;
            index = std::vector<unsigned char>();
            index.resize(size, 0);
        }//end stCellId

	/**
         * Gets the value of the bit in position i
         *
         * @param i Bit position.
         * @return the value of the bit in position i.
         */
        char getBitValue(const size_t i, const size_t DIM) { 
            const size_t j = i + (8 * index.size()) - DIM;
            const size_t p = j / 8;
            const size_t k = 7 - j % 8;
            return (index[p] & (1 << k)) > 0;     
        }//end getBitValue

        /**
         * Inverts the value of the bit in position i
         *
         * @param i Bit position.
         */
        void invertBit(const size_t i, const size_t DIM) {  
            const size_t j = i + (8 * index.size()) - DIM;
            const size_t p = j / 8;
            const size_t k = 7 - j % 8;
            index[p] ^= (1 << k);      
        }//end invertBit

	/**
         * Operator == Compare two values.
         *
         * @param cell The stCellId to be compared.
         * @return 0 if equal, <0 if this < cell, >0  if this > cell.
         */
        int operator== (const stCellId &cell) const {  
            return this->index == cell.index;     
        }//end operator ==

	void getIndex(unsigned char *index) {
	  std::copy(this->index.begin(), this->index.end(), index);
        }

        void setIndex(const unsigned char *index) {
	  std::copy(index, index + this->index.size(), this->index.begin());
        } 

        void reset() {
	  index.clear();
        }
};

#endif //__STCELLID_H


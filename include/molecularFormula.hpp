//
//  molecularFormula.hpp
//  DTarray_pro
// -----------------------------------------------------------------------------
// Copyright 2018 Aaron maurais
// -----------------------------------------------------------------------------
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// -----------------------------------------------------------------------------
//
#pragma once

#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <utility>

#include <dtarray_pro.hpp>
#include <utils.hpp>

namespace molFormula{
	
	class Species;
	class Residue;
	class Residues;
	
	typedef std::vector<std::string> HeaderType;
	typedef std::map<std::string, Species> AtomMassMapType;
	typedef std::map<std::string, int> AtomCountMapType;
	
	bool const UNICODE_AS_DEFAULT = true;
	const std::string FORMULA_RESIDUE_ORDER [] = {"C", "(13)C", "H", "D", "Br", "Cl", "N",
		"(15)N", "O", "(18)O", "P", "S", "Se"};
	const size_t FORMULA_RESIDUE_ORDER_LEN = 13;
	const std::string N_TERM_STR = "N_term";
	const std::string C_TERM_STR = "C_term";
	const std::string BAD_AMINO_ACID = "BAD_AA_IN_SEQ";
	
	std::string getFormulaFromMap(const AtomCountMapType&, bool unicode);
	//std::string symbolToUnicode(std::string);
	
	class Species{
	private:
		double mono;
		double avg;
	public:
		Species(){
			avg = 0; mono = 0;
		}
		Species(double _avg, double _mono){
			avg = _avg; mono = _mono;
		}
		~Species() {}
		
		void operator += (const Species& _add){
			avg += _add.avg; mono += _add.mono;
		}
		void operator = (const Species& _cp){
			avg = _cp.avg; mono = _cp.mono;
		}
		template<class _Tp> Species operator * (_Tp mult){
			Species ret;
			ret.avg = avg * mult;
			ret.mono = mono * mult;
			return ret;
		}
		
		double getAvg() const{
			return avg;
		}
		double getMono() const{
			return mono;
		}
	};
	
	class Residue{
	private:
		AtomCountMapType atomCountMap;
		AtomMassMapType* atomMassMap;
		
		Species masses;
		void calcMasses();
		void removeZeros();
	public:
		Residue (){
			atomMassMap = new AtomMassMapType;
		}
		Residue (AtomMassMapType* _atomMassMap){
			atomMassMap = _atomMassMap;
		}
		Residue(molFormula::AtomMassMapType* _atomMassMap,
				const molFormula::HeaderType& _header,
				const std::vector<std::string>& _elems){
			initalize(_atomMassMap, _header, _elems);
		}
		~Residue() {}
		
		//modifers
		void initalize(AtomMassMapType* , const HeaderType&, const std::vector<std::string>&);
		
		//properties
		void combineAtomCountMap(AtomCountMapType&) const;
		int getCount(std::string) const;
		std::string getFormula(bool unicode = UNICODE_AS_DEFAULT) const{
			return getFormulaFromMap(atomCountMap, unicode);
		}
		double getMass(char) const;
		double getMono() const{
			return masses.getMono();
		}
		double getAvg() const{
			return masses.getAvg();
		}
	};
	
	class Residues{
	private:
		typedef std::map<std::string, Residue> ResidueMapType;
		ResidueMapType residueMap;
		std::string atomCountTableLoc, massTableLoc;
		AtomMassMapType atomMassMap;
		HeaderType atomCountHeader;
	public:
		Residues(){
			atomCountTableLoc = ""; massTableLoc = "";
		}
		Residues(std::string _atomCountTableLoc, std::string _massTableLoc){
			atomCountTableLoc = _atomCountTableLoc; massTableLoc = _massTableLoc;
		}
		~Residues() {}
		
		bool readAtomCountTable(std::string);
		bool readAtomCountTable();
		bool readAtomMassTable(std::string);
		bool readAtomMassTable();
		bool initalize();
		bool initalize(std::string, std::string);
		
		std::string calcFormula(std::string, bool unicode = UNICODE_AS_DEFAULT,
								bool _nterm = true, bool _cterm = true) const;
		double calcMass(std::string _seq, char avg_mono,
					   bool _nterm = true, bool _cterm = true) const;
		double calcMW(std::string _seq, bool _nterm = true, bool _cterm = true) const{
			return calcMass(_seq, 'a', _nterm, _cterm);
		}
		double calcMono(std::string _seq, bool _nterm = true, bool _cterm = true) const{
			return calcMass(_seq, 'm', _nterm, _cterm);
		}
	};
}

/* molecularFormula_hpp */


//
//  BinTree.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/9/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef BinTree_hpp
#define BinTree_hpp

#include<iostream>

using namespace std;

namespace binTree{
	
	/**********************/
	/* class definitions */
	/*********************/
	
	template <class _Tp> class BinTree;
	template <class _Tp> class Node;
	
	template <class _Tp>
	class Node{
		friend class BinTree<_Tp>;
	public:
		_Tp val;
		
		Node();
		
		_Tp* getVal() const;
	private:
		Node *left;
		Node *right;
	};
	
	template <class _Tp>
	class BinTree{
	public:
		//constructor
		BinTree();
		~BinTree();
		
		//modifers
		void insert(const _Tp&);
		void destroyTree();
		
		//properties
		typename Node<_Tp>::Node* search(const _Tp&) const;
		
	private:
		typename Node<_Tp>::Node* root;
		
		//modifers
		void insert(const _Tp&, typename Node<_Tp>::Node*);
		void destroyTree(typename Node<_Tp>::Node*);
		
		//properties
		typename Node<_Tp>::Node* search(const _Tp&, typename Node<_Tp>::Node*) const;
	};
}

#endif

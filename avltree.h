/****************************************************************
Copyright (C) 2007 David M. Gay

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author disclaims all warranties with regard to this software,
including all implied warranties of merchantability and fitness.
In no event shall the author be liable for any special, indirect or
consequential damages or any damages whatsoever resulting from loss of
use, data or profits, whether in an action of contract, negligence or
other tortious action, arising out of or in connection with the use or
performance of this software.
****************************************************************/

#include <stddef.h> /* for size_t */

 typedef struct AVL_Tree AVL_Tree;
 typedef struct Element Element;

 typedef int (*AVL_Elcomp)(void*, const Element*, const Element*);
 typedef int (*AVL_Visitor)(void*, const Element*);

 extern AVL_Tree *AVL_Tree_alloc(void*, AVL_Elcomp, void *(Mallocfunc)(size_t));
 extern AVL_Tree *AVL_Tree_alloc2(void*, AVL_Elcomp, void *(Mallocfunc)(size_t), void (*Free)(void*));
 extern void AVL_Tree_free(AVL_Tree**);
 extern size_t AVL_Tree_size(AVL_Tree*);
 extern const Element *AVL_find(const Element *, AVL_Tree*);
 extern const Element *AVL_insert(const Element *, AVL_Tree*);
 extern void AVL_visit(void*, AVL_Tree*, AVL_Visitor);
 extern const Element *AVL_delete(const Element *, AVL_Tree*);
 extern void *AVL_setv(AVL_Tree *, void*);

 /* The third argument to Avl_Tree_alloc is a malloc-like function that    */
 /* only returns nonzero values.  It should use longjmp to avoid returning */
 /* if no memory is available.  If you are using the AMPL/Solver interface */
 /* library, simply pass mymalloc_ASL for this argument. */

 /* AVL_Tree_alloc returns a pointer to an AVL structure in which it has */
 /* stored its arguments.  The current void* value (possibly a "this"    */
 /* pointer) is passed to the comparision function.   AVL_setv replaces  */
 /* the current such void* value and returns the old value. */

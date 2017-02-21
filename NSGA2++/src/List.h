/*
 * List.h
 *
 *  Created on: 29/10/2009
 *      Author: antonio
 */

#ifndef LIST_H_
#define LIST_H_

typedef struct nodeList
{
    int index;
    struct nodeList *parent;
    struct nodeList *child;
}
myList;

/*
class List {
public:
   int index;
   struct lists *parent;
   struct lists *child;
};
*/

#endif /* LIST_H_ */

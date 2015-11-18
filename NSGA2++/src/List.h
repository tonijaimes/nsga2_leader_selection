/*
 * List.h
 *
 *  Created on: 29/10/2009
 *      Author: antonio
 */

#ifndef LIST_H_
#define LIST_H_

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
}
list;

/*
class List {
public:
   int index;
   struct lists *parent;
   struct lists *child;
};
*/

#endif /* LIST_H_ */

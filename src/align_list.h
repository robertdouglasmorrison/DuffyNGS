/*
 *	File		: align_list.h
 *	Content		: Double linked list which contains bam1_t align structs
 *
 * 	Created on	: 25.01.2012
 *      Author		: Wolfgang Kaisers
 *
 *	Changelog	:
 *			01.Nov.12 [get_const_next_align] Function added (returns align without copying).
 */

#ifndef ALIGN_LIST_H_
#define ALIGN_LIST_H_
#include "samtools/sam.h"
#include "samtools/bam.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// basic definitions

typedef struct align_element
{
	bam1_t *align;
	struct align_element *last_el;
	struct align_element *next_el;
} align_element;

typedef struct {
	align_element *first_el;
	align_element *last_el;
	align_element *curr_el;
	unsigned long size;
} align_list;

///////////////////////////////////////////////////////////////////////////////////////////////////
// basic functions

align_list * init_align_list()
{
	return (align_list*) calloc(1,sizeof(align_list));
}

static R_INLINE void copy_align(bam1_t *target,const bam1_t * const source)
{
	// see bam.h duplicate_align
	if(target==NULL)
		return;
	*target=*source;
	target->m_data=source->data_len;
	if ( target->data != 0) free(target->data);
	target->data=(uint8_t*)calloc(target->data_len,1);
	memcpy(target->data,source->data,target->data_len);
}

static R_INLINE bam1_t *duplicate_align(const bam1_t *src)
{
	bam1_t *b;
	b = bam_init1();
	*b = *src;
	b->m_data = b->data_len;
	b->data = (uint8_t*)calloc(b->data_len, 1);
	memcpy(b->data, src->data, b->data_len);
	return b;
}

static R_INLINE align_element* align_list_init_elem(const bam1_t *align)
{
	align_element *e=calloc(1,sizeof(align_element));
	e->align=duplicate_align(align);
	return e;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// list generic accessor functions

void align_list_push_back(align_list *l, const bam1_t *align)
{
	align_element *e=align_list_init_elem(align);
	if(l->size==0)	// list is empty
	{
		l->first_el=e;
		l->last_el=e;
		l->size=1;
	}
	else
	{
		e->last_el=l->last_el;
		e->last_el->next_el=e;
		l->last_el=e;
		++(l->size);
	}
	//printf("push_back pos: %i\n",align->core.pos);
}

void align_list_push_front(align_list *l,const bam1_t *align)
{
	align_element *e=align_list_init_elem(align);
	if(l->first_el==0)	// list is empty
	{
		l->first_el=e;
		l->last_el=e;
		l->size=1;
	}
	else
	{
		e->next_el=l->first_el;
		e->next_el->last_el=e;
		l->first_el=e;
		++(l->size);
	}
}

void align_list_pop_back(align_list *l)
{
	if(l->first_el!=l->last_el)
	{
		align_element *e=l->last_el;
		e->last_el->next_el=0;
		l->last_el=e->last_el;
		if (e) {
			if ((e->align)->data != 0) free((e->align)->data);
			free(e->align);
			free(e);
		}
		--(l->size);
	}
	else if(l->last_el>0)
	{
		//printf("pop_back pos: %i\n",l->first_el->align->core.pos);
		if (l->first_el->align->data != 0) free(l->first_el->align->data);
		free(l->first_el->align);
		free(l->first_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}

}

void align_list_pop_front(align_list *l)
{
	//printf("[pop_front] size %lu\tpos %i\n",l->size,l->first_el->align->core.pos);
	align_element *e;
	if(l->first_el!=l->last_el)
	{
		e=l->first_el;
		//printf("[pop_front] free element %lu\n",(unsigned long)e);
		if (e) {
			e->next_el->last_el=0;
			l->first_el=e->next_el;
			if ((e->align)->data != 0) free((e->align)->data);
			free(e->align);
			free(e);
		}
		--(l->size);
	}
	else if(l->first_el>0)
	{
		if (l->first_el->align->data != 0) free(l->first_el->align->data);
		free(l->first_el->align);
		free(l->first_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}
}

void wind_back(align_list *l)
{
	l->curr_el=NULL;
	return;
}


void align_list_mv_curr_elem(align_list *src,align_list *target)
{
	// Moves current element from src to end of target list
	// and moves curr_el pointer to next align

	// Nothing to do
	if(src->curr_el==NULL)
		return;

	align_element *e=src->curr_el;
	src->curr_el=e->next_el;

	///////////////////////////////////
	// Remove e from src list
	if(e->next_el!=NULL)
		e->next_el->last_el=e->last_el;
	if(e->last_el!=NULL)
		e->last_el->next_el=e->next_el;

	///////////////////////////////////
	// Insert e into end of target list
	if(target->size==0)
	{
		target->first_el=e;
		target->last_el=e;
		target->size=1;
	}
	else
	{
		target->last_el->next_el=e;
		e->last_el=target->last_el;
		target->last_el=e;
		++(target->size);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// higher level convenience functions

void destroy_align_list(align_list *l)
{
	if ( l == 0) return;
	while(l->size>0)
		align_list_pop_front(l);
	free(l);
}

bam1_t * get_next_align(align_list *l)		// Returns a *COPY* of current align
{

	if(l->first_el==NULL)
	{
		//printf("[get_next_align] l->first_el==NULL\n");
		return (bam1_t*) NULL;
	}

	if(l->curr_el==NULL)
	{
		//printf("[get_next_align] curr_el==NULL!\n");
		l->curr_el=l->first_el;
		return duplicate_align(l->curr_el->align);	// Copy!
	}
	if(l->curr_el->next_el==NULL)
	{
		//printf("[get_next_align] l->curr_el->next_el==NULL\n");
		l->curr_el=NULL;
		return (bam1_t*) NULL;
	}

	l->curr_el=l->curr_el->next_el;
	return duplicate_align(l->curr_el->align);		// Copy!
}

const bam1_t * get_const_next_align(align_list *l)		// Returns a *CONSTANT REFERENCE* to current align
{

	if(l->first_el==NULL)
	{
		return (bam1_t*) NULL;
	}

	if(l->curr_el==NULL)
	{
		l->curr_el=l->first_el;
		return l->curr_el->align;	// No Copy!
	}
	if(l->curr_el->next_el==NULL)
	{
		l->curr_el=NULL;
		return (bam1_t*) NULL;
	}

	l->curr_el=l->curr_el->next_el;
	return l->curr_el->align;		// No copy!
}


bam1_t * get_prev_align(align_list *l)
{
	if((l->last_el)==NULL)
		return (bam1_t*) NULL;
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->last_el;
		return duplicate_align(l->curr_el->align);	// Copy!
	}
	if((l->curr_el->last_el)==NULL)
	{
		l->curr_el=NULL;
		return (bam1_t*) NULL;
	}
	l->curr_el=l->curr_el->last_el;
	return duplicate_align(l->curr_el->align);
}

void write_current_align(align_list *l,bam1_t *align)
{
	if((l->curr_el)!=NULL)
		copy_align(l->curr_el->align,align);
	return;
}

void pp_curr_align(align_list *l)
{
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->first_el;
		return;
	}
	l->curr_el=(l->curr_el->next_el);
}

void mm_curr_align(align_list *l)
{
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->last_el;
		return;
	}
	l->curr_el=(l->curr_el->last_el);
}

void insert_past_curr_align(align_list *l,bam1_t *align)
{
	align_element *e=align_list_init_elem(align);
	if(l->first_el==NULL)
	{
		l->first_el=e;
		l->last_el=e;
		(l->size)=1;
		return;
	}
	if(l->curr_el==NULL)
	{
		l->first_el->last_el=e;
		e->next_el=l->first_el;
		l->first_el=e;
		++(l->size);
		return;
	}
	if(l->curr_el->next_el==NULL)
	{
		l->curr_el->next_el=e;
		e->last_el=l->curr_el;
		l->last_el=e;
		++(l->size);
		return;
	}
	e->next_el=l->curr_el->next_el;
	e->last_el=l->curr_el;
	l->curr_el->next_el->last_el=e;
	l->curr_el->next_el=e;
	++(l->size);

}

void insert_pre_curr_align(align_list *l,bam1_t *align)
{
	align_element *e=align_list_init_elem(align);
	if(l->first_el==NULL)
	{
		l->first_el=e;
		l->last_el=e;
		l->size=1;
		return;
	}
	if(l->curr_el==NULL)
	{
		l->last_el->next_el=e;
		e->last_el=l->last_el;
		l->last_el=e;
		++(l->size);
		return;
	}
	if(l->curr_el->last_el==NULL)
	{
		l->curr_el->last_el=e;
		e->next_el=l->curr_el;
		l->first_el=e;
		++(l->size);
		return;
	}
	e->last_el=l->curr_el->last_el;
	e->next_el=l->curr_el;
	l->curr_el->last_el->next_el=e;
	l->curr_el->last_el=e;
	++(l->size);
}

#endif

/* -*- coding: utf-8 -*- */
/* Copyright (C) 2000-2012 by George Williams */
/*
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.

 * The name of the author may not be used to endorse or promote products
 * derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fontforge-config.h>

#include "basics.h"

#include "ffglib.h"
#include "fontforgeui.h"
#include "multidialog.h"

/*static int MD_OK(GGadget *g, GEvent *e) {
    if ( e->type==et_controlevent && e->u.control.subtype == et_buttonactivate )
	*(int *) GDrawGetUserData(GGadgetGetWindow(g)) = true;
    return( true );
} */

static int md_e_h(GWindow gw, GEvent *event) {
    if ( event->type==et_close )
	*(int *) GDrawGetUserData(gw) = true;
    else
	return false;
    return true;
}

struct multi_expand {
    GGadgetCreateData *gcd;
    int is_row;
    int loc;
};

struct multi_postproc {
    GList_Glib *memlist;
    GList_Glib *expands;
};

static GGadgetCreateData *LayoutMultiDlgElem(MultiDlgElem *elemspec, struct multi_postproc *mpp) {
    int gcnt=1, lcnt=0, flcnt=1, g=0, l=0, fl=0, i;
    GGadgetCreateData *gcd, **flarray, *ftop;
    GTextInfo *label, *glistarray;

    if ( elemspec->question!=NULL ) {
	lcnt++;
	gcnt++;
    }
    if (elemspec->type==mde_string) {
	lcnt++;
	gcnt++;
	flcnt++;
    } else if (elemspec->type==mde_choice) {
	if (elemspec->checks) {
	    gcnt++;
	    lcnt += elemspec->answer_size;
	    gcnt += elemspec->answer_size+1;
	    flcnt += elemspec->answer_size;
	} else {
	    gcnt++;
	    flcnt++;
	}
    } else if (elemspec->type==mde_openpath || elemspec->type==mde_savepath ) {
	lcnt += 2;
	gcnt += 3;
	flcnt += 1;
    } else {
	ff_post_error(_("Unrecognized MultiAsk Element"), _("Cannot raise MultiAsk Dialog"));
	return NULL;
    }
    gcd = calloc(gcnt, sizeof(GGadgetCreateData));
    mpp->memlist = g_list_append(mpp->memlist, gcd);
    label = calloc(lcnt, sizeof(GTextInfo));
    mpp->memlist = g_list_append(mpp->memlist, label);
    flarray = calloc(flcnt, sizeof(GGadgetCreateData *));
    mpp->memlist = g_list_append(mpp->memlist, flarray);

    if (elemspec->type==mde_string) {
	label[l].text = (unichar_t *) copy(elemspec->dflt);
	label[l].text_is_1byte = true;
	gcd[g].gd.label = &label[l++];
	gcd[g].gd.flags = gg_enabled | gg_visible;
	gcd[g].data = elemspec;
	gcd[g].creator = GTextFieldCreate;
	flarray[fl++] = &gcd[g++];
    } else if (elemspec->type==mde_choice) {
	if (elemspec->checks) {
	    for ( i=0; i<elemspec->answer_size; ++i ) {
		MultiDlgAnswer *ans = &elemspec->answers[i];
		label[l].text = (unichar_t *) copy(ans->name);
		label[l].text_in_resource = true;
		label[l].text_is_1byte = true;
		gcd[g].gd.label = &label[l++];
		gcd[g].gd.flags = gg_enabled | gg_visible;
		if ( ans->is_default )
		    gcd[g].gd.flags |= gg_cb_on;
		gcd[g].data = ans;
		if ( elemspec->multiple ) {
		    gcd[g].creator = GTextFieldCreate;
		} else {
		    gcd[g].creator = GRadioCreate;
		    if ( i!=0 )
			gcd[g].gd.flags |= gg_rad_continueold;
		}
		flarray[fl++] = &gcd[g++];
	    }
	} else {
	    glistarray = calloc(elemspec->answer_size+1, sizeof(GTextInfo));
	    mpp->memlist = g_list_append(mpp->memlist, glistarray);
	    for ( i=0; i<elemspec->answer_size; ++i ) {
		MultiDlgAnswer *ans = &elemspec->answers[i];
		glistarray[i].text = (unichar_t *) copy(ans->name);
		glistarray[i].text_is_1byte = true;
		glistarray[i].selected = ans->is_default;
		glistarray[i].userdata = ans;
	    }
	    gcd[g].gd.u.list = glistarray;
	    gcd[g].gd.flags = gg_visible | gg_enabled;
	    if ( elemspec->multiple )
		gcd[g].gd.flags |= gg_list_multiplesel;
	    else
		gcd[g].gd.flags |= gg_list_exactlyone;
	    gcd[g].creator = GListCreate;
	    flarray[fl++] = &gcd[g++];
	}
    } else if (elemspec->type==mde_openpath || elemspec->type==mde_savepath ) {
	int fb=0;
        GGadgetCreateData **fbox = calloc(3, sizeof(GGadgetCreateData *));
	mpp->memlist = g_list_append(mpp->memlist, fbox);
	label[l].text = (unichar_t *) copy(elemspec->dflt);
	label[l].text_is_1byte = true;
	gcd[g].gd.label = &label[l++];
	gcd[g].gd.flags = gg_enabled | gg_visible;
	gcd[g].data = elemspec;
	gcd[g].creator = GTextFieldCreate;
	fbox[fb++] = &gcd[g++];
	label[l].text = (unichar_t *) "...";
	label[l++].text_is_1byte = true;
	gcd[g].gd.label = &label[l-1];
	gcd[g].gd.flags = gg_enabled | gg_visible;
	gcd[g].data = elemspec;
	gcd[g].creator = GButtonCreate;
	fbox[fb++] = &gcd[g++];
	gcd[g].gd.flags = gg_enabled | gg_visible;
	gcd[g].gd.u.boxelements = fbox;
	gcd[g].creator = GHBoxCreate;
	struct multi_expand *me = calloc(1, sizeof (struct multi_expand));
	mpp->expands = g_list_append(mpp->expands, me);
	me->gcd = &gcd[g];
	me->is_row = false;
	me->loc = 0;
	flarray[fl++] = &gcd[g++];
    }
    //gcd[g].gd.flags = gg_enabled | gg_visible | gg_flow_expand;
    gcd[g].gd.flags = gg_enabled | gg_visible | gg_flow_expand | gg_flow_ocenter | gg_flow_lvcenter;
    if ( !elemspec->align )
	gcd[g].gd.flags |= gg_flow_noalignlabel;
    gcd[g].gd.u.boxelements = flarray;
    gcd[g].creator = GFlowBoxCreate;
    ftop = &gcd[g++];

    if ( elemspec->question!=NULL ) {
	label[l].text = (unichar_t *) elemspec->question;
	label[l].text_is_1byte = true;
	gcd[g].gd.label = &label[l++];
	gcd[g].gd.flags = gg_enabled | gg_visible;
	gcd[g].creator = GLabelCreate;
	ftop->gd.label = (GTextInfo *) &gcd[g++];
    }
    return ftop;
}

static GGadgetCreateData *LayoutMultiDlgCategoryBody(MultiDlgCategory *catspec, struct multi_postproc *mpp) {
    int gcnt=1, flcnt=1, g=0, l=0, fl=0, i;
    GGadgetCreateData *gcd, **s1barray;

    s1barray = calloc(catspec->size+1, sizeof(GGadgetCreateData *));
    mpp->memlist = g_list_append(mpp->memlist, s1barray);
    for ( i=0; i<catspec->size; ++i ) {
	s1barray[i] = LayoutMultiDlgElem(catspec->elems+i, mpp);
	if ( s1barray[i]==NULL )
	    return false;
    }

    gcd = calloc(1, sizeof(GGadgetCreateData));
    mpp->memlist = g_list_append(mpp->memlist, gcd);
    gcd->gd.flags = gg_enabled | gg_visible | gg_s1_vert | gg_s1_flowalign | gg_s1_expand;
    gcd->gd.u.boxelements = s1barray;
    gcd->creator = GScroll1BoxCreate;

    return gcd;
}

int UI_Ask_Multi(const char *title, MultiDlgSpec *spec) {
    GRect pos, gsize;
    GWindow gw;
    struct multi_postproc mpp;
    GWindowAttrs wattrs;
    GGadgetCreateData gcd[4], boxes[3], *barray[9], *varray[3], **catgcd;
    GTextInfo label[6];
    int done = false, err = false;
    int b=0, g=0, l=0, i;

    if ( no_windowing_ui )
	return false;

    memset(&mpp,0,sizeof(mpp));
    catgcd = calloc(spec->size+1, sizeof(GGadgetCreateData *));
    mpp.memlist = g_list_append(mpp.memlist, catgcd);
    for ( i=0; i<spec->size; ++i ) {
	catgcd[i] = LayoutMultiDlgCategoryBody(spec->categories+i, &mpp);
	if ( catgcd[i]==NULL ) {
	    for ( GList_Glib *e = mpp.expands; e!=NULL; e=e->next )
		free(e->data);
	    g_list_free(mpp.expands);
	    for ( GList_Glib *m = mpp.memlist; m!=NULL; m=m->next )
		free(m->data);
	    g_list_free(mpp.memlist);
	   return false;
	}
    }

    memset(&wattrs,0,sizeof(wattrs));
    wattrs.mask = wam_events|wam_cursor|wam_utf8_wtitle|wam_undercursor|wam_isdlg|wam_restrict;
    wattrs.event_masks = ~(1<<et_charup);
    wattrs.restrict_input_to_me = 1;
    wattrs.undercursor = 1;
    wattrs.cursor = ct_pointer;
    wattrs.utf8_window_title = title;
    wattrs.is_dlg = true;
    pos.x = pos.y = 0;
    pos.width = GDrawPointsToPixels(NULL,400);
    pos.height = GDrawPointsToPixels(NULL,300);
    gw = GDrawCreateTopWindow(NULL,&pos,md_e_h,&done,&wattrs);

    memset(&label,0,sizeof(label));
    memset(&gcd,0,sizeof(gcd));
    memset(&boxes,0,sizeof(boxes));

    barray[b++] = GCD_Glue;
    label[l].text = (unichar_t *) _("_OK");
    label[l].text_is_1byte = true;
    label[l].text_in_resource = true;
    gcd[g].gd.flags = gg_visible | gg_enabled | gg_but_default;
    gcd[g].gd.label = &label[l++];
    gcd[g].creator = GButtonCreate;
    barray[b++] = &gcd[g++];
    barray[b++] = GCD_Glue;

    if ( true ) {
	label[l].text = (unichar_t *) _("_Apply");
	label[l].text_is_1byte = true;
	label[l].text_in_resource = true;
	gcd[g].gd.flags = gg_visible;
	gcd[g].gd.label = &label[l++];
	gcd[g].creator = GButtonCreate;
	barray[b++] = &gcd[g++];
	barray[b++] = GCD_Glue;
    }

    label[l].text = (unichar_t *) _("_Cancel");
    label[l].text_is_1byte = true;
    label[l].text_in_resource = true;
    gcd[g].gd.flags = gg_visible | gg_enabled | gg_but_cancel;
    gcd[g].gd.label = &label[l++];
    gcd[g].creator = GButtonCreate;
    barray[b++] = &gcd[g++];
    barray[b++] = GCD_Glue;
    barray[b++] = NULL;

    boxes[2].gd.flags = gg_enabled | gg_visible;
    boxes[2].gd.u.boxelements = barray;
    boxes[2].creator = GHBoxCreate;

    varray[0] = catgcd[0];
    varray[1] = &boxes[2];
    varray[2] = NULL;

    boxes[0].gd.flags = gg_enabled | gg_visible;
    boxes[0].gd.u.boxelements = varray;
    boxes[0].creator = GVBoxCreate;

    GGadgetsCreate(gw,boxes);

    GHVBoxSetExpandableRow(boxes[0].ret, 0);
    GHVBoxSetExpandableCol(boxes[2].ret, gb_expandglue);
    for ( GList_Glib *e = mpp.expands; e!=NULL; e=e->next ) {
	struct multi_expand *me = (struct multi_expand *) e->data;
	if ( me->is_row )
	    GHVBoxSetExpandableRow(me->gcd->ret, me->loc);
	else
	    GHVBoxSetExpandableCol(me->gcd->ret, me->loc);
	free(me);

    }
    g_list_free(mpp.expands);

    GHVBoxFitWindow(boxes[0].ret);

    for ( GList_Glib *m = mpp.memlist; m!=NULL; m=m->next )
	free(m->data);
    g_list_free(mpp.memlist);

    GDrawSetVisible(gw,true);

    while ( !done )
	GDrawProcessOneEvent(NULL);

    GDrawDestroyWindow(gw);
    return true;
}

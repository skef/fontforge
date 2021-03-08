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

static GGadgetCreateData *LayoutMultiDlgCategory(MultiDlgCategory *category, GLib_GList **memhandle) {
}

int UI_Ask_Multi(const char *title, MultiDlgSpec *spec) {
    GRect pos, gsize;
    GWindow gw;
    GLib_GList *memlist;
    GWindowAttrs wattrs;
    GGadgetCreateData gcd[12], boxes[3], *harray[3], *farray[10];
    GTextInfo label[12];
    int done = false, err = false;
    int k;

    if ( no_windowing_ui )
	return false;

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

    k = 0;
    label[k].text = (unichar_t *) _("Opt _1");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[0] = &gcd[k-1];

    label[k].text = (unichar_t *) _("Optoooooooooion _2");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[1] = &gcd[k-1];

    label[k].text = (unichar_t *) _("Optoooion _3");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[2] = &gcd[k-1];

    label[k].text = (unichar_t *) _("n _4");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[3] = &gcd[k-1];

    label[k].text = (unichar_t *) _("tion _5");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[4] = &gcd[k-1];

    label[k].text = (unichar_t *) _("xxxxxxxxOption _6");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[5] = &gcd[k-1];

    label[k].text = (unichar_t *) _("Optxxxion _7");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[6] = &gcd[k-1];

    label[k].text = (unichar_t *) _("n _8");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[7] = &gcd[k-1];

    label[k].text = (unichar_t *) _("Option _9");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GCheckBoxCreate;
    farray[8] = &gcd[k-1];
    farray[9] = NULL;

    label[k].text = (unichar_t *) _("_Options:");
    label[k].text_is_1byte = true;
    label[k].text_in_resource = true;
    gcd[k].gd.label = &label[k];
    gcd[k].gd.flags = gg_enabled | gg_visible;
    gcd[k++].creator = GLabelCreate;

    boxes[2].gd.label = (GTextInfo *) &gcd[k-1];
    boxes[2].gd.flags = gg_enabled | gg_visible;
    boxes[2].gd.u.boxelements = farray;
    boxes[2].creator = GFlowBoxCreate;
    //boxes[2].creator = GHBoxCreate;

    harray[0] = &boxes[2];
    harray[1] = NULL;

    boxes[0].gd.flags = gg_enabled | gg_visible;
    boxes[0].gd.u.boxelements = harray;
    boxes[0].creator = GScroll1BoxCreate;

    GGadgetsCreate(gw,boxes);
/*    gsize.x = gsize.y = 0;
    gsize.width = GGadgetScale(GDrawPointsToPixels(NULL,200));
    gsize.height = -1;
    GGadgetSetDesiredSize(boxes[2].ret, NULL, &gsize); */
    //GFlowBoxSetPadding(boxes[2].ret, GDrawPointsToPixels(NULL,20), -1, -1);
    GFlowBoxSetPadding(boxes[2].ret, 0, 0, -1);
    GScroll1BoxFitWindow(boxes[0].ret);

    GDrawSetVisible(gw,true);

    while ( !done )
	GDrawProcessOneEvent(NULL);

    GDrawDestroyWindow(gw);
    return true;
}

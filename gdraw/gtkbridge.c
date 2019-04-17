#include "basics.h"

#ifdef FONTFORGE_CAN_USE_GTK_BRIDGE
#include <gtk/gtk.h>
#include <gmodule.h>

struct ffgtkb_state {
	GtkWindow *event_trap_win, *the_window;
	bool grabbed, grabbing;
};
typedef struct ffgtkb_state FFGtkBState;

static void print_hello(GtkWidget *widget, gpointer data) {
	printf ("Hello World\n");
}

bool gtkb_Grabbed(FFGtkBState *state) {
	if (state == NULL)
		return false;

	return state->grabbed;
}

void gtkb_Grab(FFGtkBState *state, bool grab) {
	if (state == NULL || state->grabbing == grab)
		return;

	if (state->grabbing) {
		gtk_grab_remove(GTK_WIDGET(state->event_trap_win));
		state->grabbing = false;
		printf("Ungrabbing GTK\n");
	} else {
		gtk_grab_add(GTK_WIDGET(state->event_trap_win));
		state->grabbing = true;
		printf("Grabbing GTK\n");
	}
}

static void gtkb_GrabTrack(GtkWidget *widget, gboolean was_grabbed, gpointer data) {
	FFGtkBState *state = (FFGtkBState *) data;

	if (state == NULL) // Should never happen
		return;

	printf("notify_grab: %s\n", was_grabbed ? "True" : "False");
	state->grabbed = ! was_grabbed;
}

static void gtkb_TheWindowDestroy(GtkWidget *widget, gpointer data) {
	FFGtkBState *state = (FFGtkBState *) data;

	if (state == NULL) // Should never happen
		return;

	state->the_window = NULL;
}

FFGtkBState *gtkb_CreateState() {
	FFGtkBState *state;
	GtkWidget *window;
	GtkWidget *eb;

	state = (FFGtkBState *)calloc(1, sizeof(FFGtkBState));

	// Allocate event-trap window (will never be made visible)
	window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	eb = gtk_event_box_new();
	gtk_widget_set_sensitive(eb, true);
	gtk_container_add(GTK_CONTAINER (window), eb);
	g_signal_connect(eb, "grab-notify", G_CALLBACK (gtkb_GrabTrack), (gpointer) state);
	gtk_widget_realize(window);

	state->event_trap_win = GTK_WINDOW(window);
	return state;
}

void gtkb_DestroyState(FFGtkBState *state) {

	if (state->the_window != NULL) {
		gtk_widget_destroy(GTK_WIDGET(state->the_window));
		state->the_window = NULL;
	}

	gtk_widget_destroy(GTK_WIDGET(state->event_trap_win));
	state->event_trap_win = NULL;

	free(state);
}

static void modal_win(GtkWidget *widget, gpointer data) {
	GtkWidget *dialog, *label, *content_area;

	GtkDialogFlags flags = GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT;
	dialog = gtk_dialog_new_with_buttons("My dialog", GTK_WINDOW(gtk_widget_get_toplevel(widget)), flags,
	         "_OK", GTK_RESPONSE_ACCEPT, "_Cancel", GTK_RESPONSE_REJECT, NULL);

	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	label = gtk_label_new("This is a message");
	g_signal_connect_swapped(dialog, "response", G_CALLBACK(gtk_widget_destroy), dialog);
	gtk_container_add(GTK_CONTAINER(content_area), label);
	gtk_widget_show_all(dialog);
}

// Test window 
static GtkWidget *make_window(FFGtkBState *state, GdkWindow *w) {
	GtkWidget *window;
	GtkWidget *vbox, *view, *button, *button_box;

	window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title (GTK_WINDOW (window), "Window");
	gtk_window_set_default_size (GTK_WINDOW (window), 200, 200);
	g_signal_connect(window, "destroy", G_CALLBACK (gtkb_TheWindowDestroy), (gpointer) state);
	gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);

	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
	view = gtk_text_view_new();
	gtk_box_pack_start(GTK_BOX(vbox), view, TRUE, TRUE, 0);

	button_box = gtk_button_box_new (GTK_ORIENTATION_HORIZONTAL);
	gtk_box_pack_start(GTK_BOX(vbox), button_box, TRUE, TRUE, 0);


	button = gtk_button_new_with_label ("Open Dialog");
	//g_signal_connect (button, "clicked", G_CALLBACK (print_hello), NULL);
	g_signal_connect (button, "clicked", G_CALLBACK (modal_win), NULL);
	// g_signal_connect_swapped (button, "clicked", G_CALLBACK (gtk_widget_destroy), window);
	gtk_container_add (GTK_CONTAINER (button_box), button);

	gtk_container_add(GTK_CONTAINER (window), vbox);

	gtk_widget_show_all (window);
	// gdk_window_set_group(gtk_widget_get_window(window), gdk_window_get_group(w));
	return window;
}

void _gtkb_AddWindow(FFGtkBState *state, GdkWindow *pseudo_parent) {
	if (state != NULL && state->the_window == NULL) {
		state->the_window = GTK_WINDOW(make_window(state, pseudo_parent));
	}
}

void gtkb_ProcessEvent(FFGtkBState *state, GdkEvent *event) {
	if (state != NULL)
		gtk_main_do_event(event);
}

void gtkb_SetDefaultIcon(FFGtkBState *state, GdkPixbuf *pb) {
	if (pb != NULL) {
		GList ent = { .data = pb };
		gtk_window_set_default_icon_list(&ent);
		g_object_unref(pb);
	}
}

#endif // FONTFORGE_CAN_USE_GTK_BRIDGE


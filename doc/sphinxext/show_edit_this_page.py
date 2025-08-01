"""Sphinx extension for checking if Edit this Page button should show on this page."""

def html_page_context(app, pagename, templatename, context, doctree):
    if (
        any(part in context['exclude_edit_this_page_directory'] for part in pagename.split("/")) or
        pagename in context['exclude_edit_this_page_pagename']
    ):
        context['show_edit_this_page'] = False
    else:
        context['show_edit_this_page'] = True
            
    
def setup(app):
    app.connect('html-page-context', html_page_context)
    return {"parallel_read_safe": True}
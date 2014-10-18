
from Tkinter import *
import tkMessageBox
import ScrolledText
import sys
import tkFileDialog

my_font = ("Helvetica", "14", "bold italic")

padding = 15
ipadding = 2

def run_analysis(c1, c2, a, l, f):
    print "Fibration choice: %s" % f
    print "Analysis type: %s" % a
    print "Logging level: %s" % l
    print "Showing graphics: %s" % c1
    print "Saving to file: %s" % c2

class RedirectText(object):
    """"""

    def __init__(self, text_ctrl):
        """Constructor"""
        self.output = text_ctrl

    def write(self, string):
        """"""
        self.output.insert(END, string)


class Application(Frame):

    def main_interface(self):
        new_window = Toplevel()
        new_window.geometry("800x600")
        new_window.wm_title("New Analysis")

        analysis_var = StringVar()
        analysis_var.set("no analysis chosen")
        single_network_button = Radiobutton(new_window,
                                text="single network", 
                                font=my_font,
                                variable = analysis_var,
                                value = "single"
                                )
        phase_scan_button = Radiobutton(new_window,
                                text="phase scan", 
                                font=my_font,
                                variable = analysis_var,
                                value = "full"
                                )

        log_var = StringVar()
        log_var.set("no log_level chosen")
        log_warning_button = Radiobutton(new_window,
                                text="warning-level logging", 
                                font=my_font,
                                variable = log_var,
                                value = "warning"
                                )
        log_info_button = Radiobutton(new_window,
                                text="info-level logging", 
                                font=my_font,
                                variable = log_var,
                                value = "info"
                                )
        log_debug_button = Radiobutton(new_window,
                                text="debug-level logging", 
                                font=my_font,
                                variable = log_var,
                                value = "debug"
                                )

        g_var = BooleanVar()
        g_var.set(False)
        g_button = Checkbutton(new_window, 
                        text = "Show Graphics", 
                        variable = g_var, \
                        onvalue = True, 
                        offvalue = False, 
                        # height=2, \
                        # width = 20,
                        font = my_font,
                        )

        w_var = BooleanVar()
        w_var.set(False)
        w_button = Checkbutton(new_window, 
                        text = "Save to file", 
                        variable = w_var, \
                        onvalue = True, 
                        offvalue = False, 
                        # height=2, \
                        # width = 20,
                        font = my_font
                        )

        run = Button(new_window,
                        text = "Run!",
                        height = 2,
                        width = 20,
                        command = lambda: \
                            run_analysis( \
                                g_var.get(), \
                                w_var.get(), \
                                analysis_var.get(), \
                                log_var.get(), \
                                fibration_list.get(\
                                    fibration_list.curselection() \
                                    ) \
                                ),
                        font = my_font
                        )

        close = Button(new_window,
                        text = "Close",
                        height = 2,
                        width = 20,
                        command = new_window.destroy,
                        font = my_font
                        )

        terminal_frame = Frame(new_window)
        terminal_output = ScrolledText.ScrolledText(
                                                    terminal_frame, 
                                                    height = 40,
                                                    width = 80
                                                    )
        # redirect stdout
        redir = RedirectText(terminal_output)
        sys.stdout = redir

        fibration_list = self.fetch_fibrations(new_window)
        fibration_label = Label(new_window, text="Fibrations", font = my_font)

        separator_0 = Frame(new_window, height=2, bd=1, relief=SUNKEN)
        separator_1 = Frame(new_window, height=2, bd=1, relief=SUNKEN)
        separator_2 = Frame(new_window, height=2, bd=1, relief=SUNKEN)
        separator_3 = Frame(new_window, height=2, bd=1, relief=SUNKEN)
        separator_4 = Frame(new_window, height=2, bd=1, relief=SUNKEN)

        # defaults
        single_network_button.select()
        log_info_button.select()
        g_button.select()
        w_button.deselect()
        fibration_list.select_set(0)

        # organize
        terminal_frame.pack(anchor='ne', side="right", fill=Y)
        terminal_output.pack()
        # separator.pack(padx=5, pady=5, anchor='ne', side=right)
        #
        fibration_label.pack(anchor='nw',side="top")
        fibration_list.pack(anchor='nw',side="top")
        #
        separator_0.pack(fill=X, padx=5, pady=5, anchor='nw', side="top")
        #
        single_network_button.pack(anchor = 'nw',side="top")
        phase_scan_button.pack(anchor = 'nw',side="top")
        #
        separator_1.pack(fill=X, padx=5, pady=5, anchor='nw', side="top")
        #
        g_button.pack(anchor = 'nw',side="top")
        w_button.pack(anchor = 'nw',side="top")
        #
        separator_2.pack(fill=X, padx=5, pady=5, anchor='nw', side="top")
        #
        log_warning_button.pack(anchor = 'nw',side="top")
        log_info_button.pack(anchor = 'nw',side="top")
        log_debug_button.pack(anchor = 'nw',side="top")
        #
        separator_3.pack(fill=X, padx=5, pady=15, anchor='nw', side="top")
        #
        run.pack(anchor = 'nw',side="top")
        close.pack(anchor = 'nw',side="top")
        


    def fetch_fibrations(self, window):
        lb = Listbox(window)
        lb.insert(1, "su2")
        lb.insert(2, "invented")
        lb.insert(3, "Nf=1")
        lb.insert(4, "another one!")
        return lb

    def load_old(self):
        tkMessageBox.showinfo( "Work In Progress", 
                                "This feature is not yet available!")

    def createWidgets(self):
        self.QUIT = Button(self, font=my_font)
        self.QUIT["text"] = "QUIT"
        self.QUIT["fg"]   = "red"
        self.QUIT["command"] =  self.quit

        self.QUIT.pack(side="bottom", 
                        padx = padding, 
                        pady = padding,
                        ipadx = ipadding,
                        ipady = ipadding
                        )

        self.new_analysis = Button(self, text="New analysis", font=my_font)
        self.new_analysis["command"] = self.main_interface
        self.new_analysis.pack(side="top", 
                        padx = padding, 
                        pady = padding,
                        ipadx = ipadding,
                        ipady = ipadding
                        )

        self.load_analysis = Button(self, text="Load analysis", font=my_font)
        self.load_analysis["command"] = self.load_old
        self.load_analysis.pack(side="bottom", 
                        padx = padding, 
                        pady = padding,
                        ipadx =ipadding,
                        ipady = ipadding
                        )


    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

root = Tk()
root.geometry("300x200")
root.wm_title("MOSE - MOduli Space Explorer")
app = Application(master=root)
app.mainloop()
root.destroy()
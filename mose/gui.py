from Tkinter import *
import tkMessageBox
import ScrolledText
import sys
import tkFileDialog
import os

# import tkSimpleDialog

from api import analysis

my_font = ("Helvetica", "14", "bold italic")

padding = 15
ipadding = 2

def run_analysis(window, graphics, save, analysis_type, log, fibration):
    print "Fibration choice: %s" % fibration
    print "Analysis type: %s" % analysis_type
    print "Logging level: %s" % log
    print "Showing graphics: %s" % str(graphics)
    print "Saving to file: %s" % str(save)
    if analysis_type == 'single':
        ask_phase = PhaseDiag(window)
        window.wait_window(ask_phase.top)
        phase = ask_phase.var
        analysis(graphics, save, analysis_type, log, fibration, phase=phase)
    else:
        ask_range = PhaseRangeDiag(window)
        window.wait_window(ask_range.top)
        theta_range = ask_range.var
        analysis(graphics, save, analysis_type, log, fibration, 
                 theta_range=theta_range,)

class STDText(Text):
    def __init__(self, parent):
        Text.__init__(self, parent)
        self.parent=parent

    def write(self, stuff):
        self.insert("end", stuff)
        self.yview_pickplace("end")

    def flush(self):
        None
        


class PhaseDiag:

    def __init__(self, parent):

        top = self.top = Toplevel(parent)

        Label(top, text="Phase of the K-wall Network:").pack()
        self.e = Entry(top)
        self.e.pack(padx=5)
        self.var = None

        b = Button(top, text="OK", command=self.ok)
        b.pack(pady=5)

    def ok(self):

        # print "value is", self.e.get()
        self.var = float(self.e.get())

        self.top.destroy()


class PhaseRangeDiag:
    def __init__(self, parent):

        top = self.top = Toplevel(parent)

        
        self.theta_0 = Entry(top)
        self.theta_1 = Entry(top)
        self.steps = Entry(top)
        Label(top, text="Initial phase").pack(side="top")
        self.theta_0.pack(padx=5, side='top')
        Label(top, text="Final phase").pack(side="top")
        self.theta_1.pack(padx=5,  side='top')
        Label(top, text="Steps").pack(side="top")
        self.steps.pack(padx=5, side='top')
        self.var = None

        b = Button(top, text="OK", command=self.ok)
        b.pack(pady=5)

    def ok(self):

        # print "value is", self.e.get()
        self.var = [
                    float(self.theta_0.get()),
                    float(self.theta_1.get()),
                    int(self.steps.get())
                    ]

        self.top.destroy()


def walk_dir(root_dir, extension):
    """
    retrieve files with specified extension 
    within the specified directory
    """
    file_list = []
    full_paths = []
    for path in os.listdir(root_dir):
            path = os.path.join(root_dir, path).lower()
            # if os.path.isfile(path) and path.endswith(extension):
            if path.endswith(extension):
                file_list.append(os.path.split(path)[1])
                # I'm not using this, but it may be useful in the future
                full_paths.append(path) 
    print "This is the folder %s" % root_dir
    if len(file_list) == 0:
        print "Did not find any fibration files in folder\n%s !" % root_dir
    return file_list


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
                        text = "show graphics", 
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
                        text = "save to file", 
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
                                new_window, \
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

        # terminal_frame = Frame(new_window)
        # terminal_output = ScrolledText.ScrolledText(
        #                                             terminal_frame, 
        #                                             height = 40,
        #                                             width = 80
        #                                             )
        # # redirect stdout
        # redir = RedirectText(terminal_output)
        # sys.stdout = redir

        terminal_frame = STDText(new_window)
        sys.stdout = terminal_frame
        sys.stderr = terminal_frame

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
        # terminal_output.pack()
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
        root_dir = os.getcwd()
        extension = '.ini'
        file_list = walk_dir(root_dir, extension)
        file_list.sort()
        for i, f_name in list(enumerate(file_list)):
            lb.insert(i, f_name)
        return lb

    def load_old(self):
        """
        load a .mose file and allow to analyze the data, 
        whatever that means..
        """
        tkMessageBox.showinfo( "Work In Progress", 
                                "This feature is not yet available!")

    def fibration_creator(self):
        """
        a mask that allows to enter parameters that will 
        be translated into a formatted .ini file, which 
        will be stored and made available for future use
        """
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

        self.create_fibration = Button( 
                                        self, 
                                        text="Create new fibration", 
                                        font=my_font
                                        )
        self.create_fibration["command"] = self.fibration_creator
        self.create_fibration.pack(side="bottom", 
                        padx = padding, 
                        pady = padding,
                        ipadx =ipadding,
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

def open_gui(config, data):
    root = Tk()
    root.geometry("300x300")
    root.wm_title("MOSE - MOduli Space Explorer")
    app = Application(master=root)
    app.mainloop()
    root.destroy()

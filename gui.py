
from Tkinter import *
import tkFont
import tkMessageBox

# helv36 = tkFont.Font(family='Helvetica', size=36, weight='bold') 
my_font = ("Helvetica", "14", "bold italic")

padding = 15
ipadding = 2

def run_analysis(c1,c2):
        print "button 1: %s" % c1
        print "button 2: %s" % c2

class Application(Frame):
    def main_interface(self):
        new_window = Toplevel()
        new_window.geometry("600x400")
        new_window.wm_title("New Analysis")


        g_var = IntVar()
        w_var = IntVar()
        g_button = Checkbutton(new_window, 
                        text = "Show Graphics", 
                        variable = g_var, \
                        onvalue = 1, 
                        offvalue = 0, 
                        height=2, \
                        width = 20,
                        font = my_font,
                        )
        w_button = Checkbutton(new_window, 
                        text = "Save to file", 
                        variable = w_var, \
                        onvalue = 1, 
                        offvalue = 0, 
                        height=2, \
                        width = 20,
                        font = my_font
                        )
        run = Button(new_window,
                        text = "Run!",
                        height = 2,
                        width = 20,
                        command = lambda: \
                            run_analysis( \
                                g_var.get(), \
                                w_var.get() \
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
        g_button.pack(anchor='w', side="top")
        w_button.pack(anchor='w', side="top")
        run.pack(anchor='w', side="top")
        close.pack(anchor='w', side="top")



    def load_old(self):
        tkMessageBox.showinfo( "Work In Progress", 
                                "This functionality is not ready yet!")

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
                        padx = 10, 
                        pady = 20,
                        ipadx =2,
                        ipady = 10
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
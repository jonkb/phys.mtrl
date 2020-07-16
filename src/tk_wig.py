import math
import tkinter as tk
from tkinter import font

#Expansion of tk.Text, used for displaying text, not receiving any input
#You can't copy text from a Label, which is a problem
class Txt_wig(tk.Text):
	ft_fam = "Helvetica"#"MS Sans Serif"#"System"
	ft_size = 11
	padw = 4
	#rep_lbl = tk.Label(root, text=report, justify=tk.LEFT)
	#tk.Text allows for copy / paste
	def __init__(self, root, txt):
		numlines = txt.count("\n") + 1
		#TO DO: One too few when there's a blank line ?????
		#print(374, numlines)
		super().__init__(root, relief="flat", height=numlines)
		
		#print(555, txt)
		self.txt = txt
		self.insert(tk.END, txt)
		#int(self.index("end-1c").split('.')[0])
		self.TWfont = font.Font(root=root, family=self.ft_fam, size=self.ft_size)
		#print(self.TWfont.actual())
		#self.TWfont = font.Font(family="System", size=12)
		#print(font.families(root=root))
		#print(self.TWfont)
		self.config(font=self.TWfont)
		#self.config(font=("System", 12))
		self.config(state="disabled")
		self.bind("<Configure>", self.resize)
	def packslf(self):
		self.pack(fill="both", padx=1, pady=1)
	def resize(self, event):
		#self.config(font=self.TWfont)
		#print(381, event.width, event.height)
		textw = event.width - self.padw
		#print(382, len(self.get("1.0", "end-1c").split("\n")))
		#print(390, self.TWfont.measure("abcdefghi"))
		#print(390.5, self.TWfont.measure("abc def ghi"))
		#print(391, self.TWfont.measure("abcde\tfghi"))
		#print(391, self.TWfont.measure("\tabcdefghi"))
		#print(391, self.TWfont.measure("\n\tabcdefghi"))
		#print(382, self.TWfont.measure(self.txt.split('\n')[0]))
		height = 0
		for line in self.txt.split('\n'):
			#Using max(1, l) means that empty lines still add to the height
			height += max(1, math.ceil(self.TWfont.measure(line) / textw))
		self.config(height=height)
	
def txtwigtest():
	root = tk.Tk()
	tw = Txt_wig(root, "ABCDEFGHIJKLMNOPQRSTUVWXYZ\nDEF\nGHIJKLM")
	tw.packslf()
	root.mainloop()
#txtwigtest()
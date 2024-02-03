import datetime as dati
import tkinter as tk
import tkinter.ttk as ttk

# main window
window = tk.Tk()
# window title
window.title("feline")
# Set a minimum window size (width x height)
window.minsize(1920, 1080)
# open window maximized
window.attributes("-zoomed", True)

label_text = tk.StringVar()
label_text.set(f"Uhrzeit: {dati.datetime.now().strftime('%H:%M:%S')}")
label = ttk.Label(window, textvariable=label_text, font=("Helvetica", 14))
button = ttk.Button(window, text="Aktuelle Zeit anzeigen", command=lambda: label_text.set(f"Uhrzeit: {dati.datetime.now().strftime('%H:%M:%S')}"))

# Center the label using pack geometry manager
label.pack(expand=True)
button.pack(expand=True)

# build window
window.mainloop()

import tkinter as tk
import datetime as dati

# main window
window = tk.Tk()
# window title
window.title("feline")
# open window maximised
window.attributes("-zoomed", True)

label_text = tk.StringVar()
label_text.set(f"Uhrzeit: {dati.datetime.now().strftime('%H:%M:%S')}")
label = tk.Label(window, textvariable=label_text, font=("Helvetica", 14))
button = tk.Button(window, text="Aktuelle Zeit anzeigen", command=lambda: label_text.set(f"Uhrzeit: {dati.datetime.now().strftime('%H:%M:%S')}"))

# Center the label using pack geometry manager
label.pack(expand=True)
button.pack()

# build window
window.mainloop()

from glidesign.geometry import LiftingSurface
    
if __name__ == '__main__':
    from parapy.gui import display
    wing = LiftingSurface(name = "Wing")
    display(wing)

class LabeledRectangle:
    lock = None  # only one can be animated at a time

    def __init__(self, rect, channel, color):
        self.rect = rect
        self.press = None
        self.background = None
        self.channel_str = str(channel)
        axes = self.rect.axes
        x0, y0 = self.rect.xy
        self.text = axes.text(x0, y0, self.channel_str, color=color)
        self.text.set_visible(False)

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes:
            return
        if LabeledRectangle.lock is not None:
            return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        LabeledRectangle.lock = self
        self.text.set_visible(True)
        self.text.draw()

    def on_release(self, event):
        'on release we reset the press data'
        if LabeledRectangle.lock is not self:
            return
        self.press = None
        LabeledRectangle.lock = None
        self.text.set_visible(False)
        self.text.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)


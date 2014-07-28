import wtforms as wtf
from parampool.html5.flask.fields import FloatField

class Main(wtf.Form):
    I            = FloatField(default=1.0,
                             validators=[wtf.validators.InputRequired()])
    a            = FloatField(default=0.2,
                             validators=[wtf.validators.InputRequired()])
    T            = FloatField(default=4.0,
                             validators=[wtf.validators.InputRequired()])
    dt_values    = wtf.TextField(default='[1.25, 0.75, 0.5, 0.1]',
                             validators=[wtf.validators.InputRequired()])
    theta_values = wtf.TextField(default='[0, 0.5, 1]',
                             validators=[wtf.validators.InputRequired()])

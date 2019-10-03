from django.contrib import admin
from keckdatabase.models import Object, Observations, Measurements

# Register your models here.
admin.site.register(Object)
admin.site.register(Observations)
admin.site.register(Measurements)
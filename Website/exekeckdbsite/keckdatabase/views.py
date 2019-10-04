from django.views import generic
from django.shortcuts import render
from keckdatabase.models import Object, Measurements, Observations

def index(request):
    num_objects = Object.objects.all().count()
    object_list = Object.objects.all()

    context = {
        'num_objects' : num_objects,
    }

    return render(request, 'index.html', context=context)

class ObjectListView(generic.ListView):
    model = Object

    def get_context_data(self, **kwargs):
        context = super(ObjectListView, self).get_context_data(**kwargs)
        context['measurement'] = Measurements.objects.all()
        context['observation'] = Observations.objects.all()

        return context

class ObjectDetailView(generic.DetailView):
    model = Object
    context_object_name = "object-detail"

    def get_context_data(self, **kwargs):
        context = super(ObjectDetailView, self).get_context_data(**kwargs)
        ob_id = self.kwargs['pk']
        context['measurement'] = Measurements.objects.get(mid=ob_id)
        context['observation'] = Observations.objects.get(obid=ob_id)

        return context
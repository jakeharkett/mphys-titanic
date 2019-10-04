from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('objects/', views.ObjectListView.as_view(), name="objects"),
    path('objects/<int:pk>', views.ObjectDetailView.as_view(), name="object-detail"),
]
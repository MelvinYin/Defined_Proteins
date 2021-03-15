from django.urls import path
from query_main import views
from django.conf.urls.static import static


urlpatterns = [path('', views.index, name='index'),
               path('results1', views.results1, name='results1'),
               path('results2', views.results2, name='results2'),
               path('results3', views.results3, name='results3'),
               ] + static('/static/tmp')
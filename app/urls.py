from django.urls import path
from app.views import (
    index,
    api
)
urlpatterns = [
    path("", index, name="index"),
    path("api/", api.urls),
]

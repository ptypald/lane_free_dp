
var app = {
    el: "#app",
    data: {
        toolbar: {
            title: "Challenge 1",
            color: "",
            dark: false
        },
		  projects: []
    }
};

Vue.use(Vuetify);
new Vue(app);

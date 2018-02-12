(function() 
{
    document.addEventListener("DOMContentLoaded", function(event) 
    {
        // Don't do anything if we're running on localhost or viewing as a file
        // The 'nbviewer' service won't work, and at least we'll be able to
        // verify that the link works
        if (location.hostname == "" || location.hostname == "localhost")
            return;

        var anchors = document.getElementsByTagName("a");
    
        for (var i = 0 ; i < anchors.length ; i++)
        {
            var a = anchors[i];
            var href = a.getAttribute("href");

            // Change links that point to notebooks so that they can be viewed
            // automatically with http://nbviewer.jupyter.org
            if (href.startsWith("_static/") && href.endsWith(".ipynb"))
            {
                console.log(href);

                var curPath = location.pathname;
                var idx = curPath.length - 1;
                while (idx >= 0 && curPath[idx] != '/')
                    idx--;

                if (idx < 0)
                    curPath = "/";
                else
                    curPath = curPath.substr(0, idx+1);

                var nbUrl = "https://nbviewer.jupyter.org/url/" + location.host + curPath + href;

                a.setAttribute("href", nbUrl);
                a.setAttribute("target", "_blank"); // Open a new window/tab
            }
        }
    });
})();

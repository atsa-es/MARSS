MathJax = {
  loader: {load: ['[tex]/tagformat', '[tex]/ams']},
  section: 1,
  tex: {
    tags: 'ams',
    packages: {'[+]': ['tagformat', 'sections', 'ams']},
    tagformat: {
      number: (n) => MathJax.config.section + '.' + n
    }
  },
  startup: {
    ready() {
      const Configuration = MathJax._.input.tex.Configuration.Configuration;
      const CommandMap = MathJax._.input.tex.SymbolMap.CommandMap;
      new CommandMap('sections', {
        nextSection: 'NextSection'
      }, {
        NextSection(parser, name) {
          MathJax.config.section++;
          parser.tags.counter = parser.tags.allCounter = 0;
        }
      });
      Configuration.create(
        'sections', {handler: {macro: ['sections']}}
      );
      MathJax.startup.defaultReady();
    }
  }
};

  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml-full.js";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();

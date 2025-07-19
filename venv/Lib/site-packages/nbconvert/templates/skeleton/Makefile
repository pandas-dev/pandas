TPLS := $(patsubst %.tpl,../latex/skeleton/%.tplx,$(wildcard *.tpl))

all: clean $(TPLS)

# Convert standard Jinja2 syntax to LaTeX safe Jinja2
# see http://flask.pocoo.org/snippets/55/ for more info
../latex/skeleton/%.tplx: %.tpl
	@echo 'generating tex equivalent of $^: $@'
	@echo '((=- Auto-generated template file, DO NOT edit directly!\n' \
		  '   To edit this file, please refer to ../../skeleton/README.md' \
		  '-=))\n\n' > $@
	@sed \
		-e 's/{%/((*/g' \
		-e 's/%}/*))/g' \
		-e 's/{{/(((/g' \
		-e 's/}}/)))/g' \
		-e 's/{#/((=/g' \
		-e 's/#}/=))/g' \
		-e "s/tpl'/tplx'/g" \
		$^ >> $@

clean:
	@echo "cleaning generated tplx files..."
	@-rm ../latex/skeleton/*.tplx

(function() {
	/* Reformats last part of 'whatsnew' URLs */
	const generateWhatsNew = (whatsNew) => {
		let url = whatsNew[0], anchor = whatsNew[1];
		url = ((url) => url.split(".")[0] + "/")(url);
		let anchorEls = anchor.slice(9).split("-");
		anchorEls[0] = ((version) => "v" + version.slice(0,1) + "." + version.slice(1,-1) + "." + version.slice(-1) + ".html#")(anchorEls[0]);
		return url + anchorEls[0] + anchorEls.slice(1).join("-");
	}

	const links = Array.from(document.getElementsByTagName("a"));
	links.forEach((link) => {
		const re = /(whatsnew.html)#(whatsnew)-[\d]+(-\w+)+/g;
		let linkElements = link.href.split("/");		
		if (re.test(linkElements.slice(-1)[0])) {
			let whatsNew = linkElements.slice(-1)[0].split("#");
			whatsNew = generateWhatsNew(whatsNew);
			linkElements[linkElements.length - 1] = whatsNew;
			link.href = linkElements.join("/")
		}
	});
})();

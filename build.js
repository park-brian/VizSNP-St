import { build } from "esbuild";
import { polyfillNode } from "esbuild-plugin-polyfill-node";

build({
	entryPoints: ["src/app.js"],
	bundle: true,
	minify: true,
	sourcemap: true,
	outfile: "docs/app.js",
	plugins: [
		polyfillNode({
      polyfills: ['buffer']
		}),
	],
});
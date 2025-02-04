const express = require("express");
const path = require("path");

const app = express();
app.use(express.json());
app.use(express.static(path.join(__dirname, "..", "public")));

app.use("/rust-ev-calc/pkg", express.static(path.join(__dirname, "..", "rust-ev-calc", "pkg")));

app.listen(3000, () => {
  console.log("Server running at http://localhost:3000");
});

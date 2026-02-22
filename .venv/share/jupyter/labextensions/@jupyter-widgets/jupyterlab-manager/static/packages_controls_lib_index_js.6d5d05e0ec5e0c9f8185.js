"use strict";
(self["webpackChunk_jupyter_widgets_jupyterlab_manager"] = self["webpackChunk_jupyter_widgets_jupyterlab_manager"] || []).push([["packages_controls_lib_index_js"],{

/***/ "../../packages/controls/lib/index.js":
/*!********************************************!*\
  !*** ../../packages/controls/lib/index.js ***!
  \********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   AccordionModel: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.AccordionModel),
/* harmony export */   AccordionView: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.AccordionView),
/* harmony export */   AudioModel: () => (/* reexport safe */ _widget_audio__WEBPACK_IMPORTED_MODULE_8__.AudioModel),
/* harmony export */   AudioView: () => (/* reexport safe */ _widget_audio__WEBPACK_IMPORTED_MODULE_8__.AudioView),
/* harmony export */   BaseIntSliderView: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.BaseIntSliderView),
/* harmony export */   BoolModel: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.BoolModel),
/* harmony export */   BoundedFloatModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.BoundedFloatModel),
/* harmony export */   BoundedFloatTextModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.BoundedFloatTextModel),
/* harmony export */   BoundedIntModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.BoundedIntModel),
/* harmony export */   BoundedIntTextModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.BoundedIntTextModel),
/* harmony export */   BoxModel: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.BoxModel),
/* harmony export */   BoxView: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.BoxView),
/* harmony export */   ButtonModel: () => (/* reexport safe */ _widget_button__WEBPACK_IMPORTED_MODULE_4__.ButtonModel),
/* harmony export */   ButtonStyleModel: () => (/* reexport safe */ _widget_button__WEBPACK_IMPORTED_MODULE_4__.ButtonStyleModel),
/* harmony export */   ButtonView: () => (/* reexport safe */ _widget_button__WEBPACK_IMPORTED_MODULE_4__.ButtonView),
/* harmony export */   CheckboxModel: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.CheckboxModel),
/* harmony export */   CheckboxStyleModel: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.CheckboxStyleModel),
/* harmony export */   CheckboxView: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.CheckboxView),
/* harmony export */   ColorPickerModel: () => (/* reexport safe */ _widget_color__WEBPACK_IMPORTED_MODULE_9__.ColorPickerModel),
/* harmony export */   ColorPickerView: () => (/* reexport safe */ _widget_color__WEBPACK_IMPORTED_MODULE_9__.ColorPickerView),
/* harmony export */   ColorsInputModel: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.ColorsInputModel),
/* harmony export */   ColorsInputView: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.ColorsInputView),
/* harmony export */   ComboboxModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.ComboboxModel),
/* harmony export */   ComboboxView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.ComboboxView),
/* harmony export */   ControllerAxisModel: () => (/* reexport safe */ _widget_controller__WEBPACK_IMPORTED_MODULE_15__.ControllerAxisModel),
/* harmony export */   ControllerAxisView: () => (/* reexport safe */ _widget_controller__WEBPACK_IMPORTED_MODULE_15__.ControllerAxisView),
/* harmony export */   ControllerButtonModel: () => (/* reexport safe */ _widget_controller__WEBPACK_IMPORTED_MODULE_15__.ControllerButtonModel),
/* harmony export */   ControllerButtonView: () => (/* reexport safe */ _widget_controller__WEBPACK_IMPORTED_MODULE_15__.ControllerButtonView),
/* harmony export */   ControllerModel: () => (/* reexport safe */ _widget_controller__WEBPACK_IMPORTED_MODULE_15__.ControllerModel),
/* harmony export */   ControllerView: () => (/* reexport safe */ _widget_controller__WEBPACK_IMPORTED_MODULE_15__.ControllerView),
/* harmony export */   DatePickerModel: () => (/* reexport safe */ _widget_date__WEBPACK_IMPORTED_MODULE_10__.DatePickerModel),
/* harmony export */   DatePickerView: () => (/* reexport safe */ _widget_date__WEBPACK_IMPORTED_MODULE_10__.DatePickerView),
/* harmony export */   DatetimeModel: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.DatetimeModel),
/* harmony export */   DatetimeView: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.DatetimeView),
/* harmony export */   DescriptionModel: () => (/* reexport safe */ _widget_description__WEBPACK_IMPORTED_MODULE_20__.DescriptionModel),
/* harmony export */   DescriptionStyleModel: () => (/* reexport safe */ _widget_description__WEBPACK_IMPORTED_MODULE_20__.DescriptionStyleModel),
/* harmony export */   DescriptionView: () => (/* reexport safe */ _widget_description__WEBPACK_IMPORTED_MODULE_20__.DescriptionView),
/* harmony export */   DirectionalLinkModel: () => (/* reexport safe */ _widget_link__WEBPACK_IMPORTED_MODULE_2__.DirectionalLinkModel),
/* harmony export */   DropdownModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.DropdownModel),
/* harmony export */   DropdownView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.DropdownView),
/* harmony export */   FileUploadModel: () => (/* reexport safe */ _widget_upload__WEBPACK_IMPORTED_MODULE_21__.FileUploadModel),
/* harmony export */   FileUploadView: () => (/* reexport safe */ _widget_upload__WEBPACK_IMPORTED_MODULE_21__.FileUploadView),
/* harmony export */   FloatLogSliderModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatLogSliderModel),
/* harmony export */   FloatLogSliderView: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatLogSliderView),
/* harmony export */   FloatModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatModel),
/* harmony export */   FloatProgressModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatProgressModel),
/* harmony export */   FloatRangeSliderModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatRangeSliderModel),
/* harmony export */   FloatRangeSliderView: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatRangeSliderView),
/* harmony export */   FloatSliderModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatSliderModel),
/* harmony export */   FloatSliderView: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatSliderView),
/* harmony export */   FloatTextModel: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatTextModel),
/* harmony export */   FloatTextView: () => (/* reexport safe */ _widget_float__WEBPACK_IMPORTED_MODULE_14__.FloatTextView),
/* harmony export */   FloatsInputModel: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.FloatsInputModel),
/* harmony export */   FloatsInputView: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.FloatsInputView),
/* harmony export */   GridBoxModel: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.GridBoxModel),
/* harmony export */   GridBoxView: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.GridBoxView),
/* harmony export */   HBoxModel: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.HBoxModel),
/* harmony export */   HBoxView: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.HBoxView),
/* harmony export */   HTMLMathModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.HTMLMathModel),
/* harmony export */   HTMLMathStyleModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.HTMLMathStyleModel),
/* harmony export */   HTMLMathView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.HTMLMathView),
/* harmony export */   HTMLModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.HTMLModel),
/* harmony export */   HTMLStyleModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.HTMLStyleModel),
/* harmony export */   HTMLView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.HTMLView),
/* harmony export */   ImageModel: () => (/* reexport safe */ _widget_image__WEBPACK_IMPORTED_MODULE_6__.ImageModel),
/* harmony export */   ImageView: () => (/* reexport safe */ _widget_image__WEBPACK_IMPORTED_MODULE_6__.ImageView),
/* harmony export */   IntModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntModel),
/* harmony export */   IntProgressModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntProgressModel),
/* harmony export */   IntRangeSliderModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntRangeSliderModel),
/* harmony export */   IntRangeSliderView: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntRangeSliderView),
/* harmony export */   IntSliderModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntSliderModel),
/* harmony export */   IntSliderView: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntSliderView),
/* harmony export */   IntTextModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntTextModel),
/* harmony export */   IntTextView: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.IntTextView),
/* harmony export */   IntsInputModel: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.IntsInputModel),
/* harmony export */   IntsInputView: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.IntsInputView),
/* harmony export */   JUPYTER_CONTROLS_VERSION: () => (/* reexport safe */ _version__WEBPACK_IMPORTED_MODULE_1__.JUPYTER_CONTROLS_VERSION),
/* harmony export */   JupyterLuminoAccordionWidget: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.JupyterLuminoAccordionWidget),
/* harmony export */   JupyterLuminoTabPanelWidget: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.JupyterLuminoTabPanelWidget),
/* harmony export */   LabelModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.LabelModel),
/* harmony export */   LabelStyleModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.LabelStyleModel),
/* harmony export */   LabelView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.LabelView),
/* harmony export */   LabeledDOMWidgetModel: () => (/* reexport safe */ _widget_description__WEBPACK_IMPORTED_MODULE_20__.LabeledDOMWidgetModel),
/* harmony export */   LabeledDOMWidgetView: () => (/* reexport safe */ _widget_description__WEBPACK_IMPORTED_MODULE_20__.LabeledDOMWidgetView),
/* harmony export */   LinkModel: () => (/* reexport safe */ _widget_link__WEBPACK_IMPORTED_MODULE_2__.LinkModel),
/* harmony export */   MultipleSelectionModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.MultipleSelectionModel),
/* harmony export */   NaiveDatetimeModel: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.NaiveDatetimeModel),
/* harmony export */   PasswordModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.PasswordModel),
/* harmony export */   PasswordView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.PasswordView),
/* harmony export */   PlayModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.PlayModel),
/* harmony export */   PlayView: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.PlayView),
/* harmony export */   ProgressStyleModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.ProgressStyleModel),
/* harmony export */   ProgressView: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.ProgressView),
/* harmony export */   RadioButtonsModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.RadioButtonsModel),
/* harmony export */   RadioButtonsView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.RadioButtonsView),
/* harmony export */   SelectModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectModel),
/* harmony export */   SelectMultipleModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectMultipleModel),
/* harmony export */   SelectMultipleView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectMultipleView),
/* harmony export */   SelectView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectView),
/* harmony export */   SelectionContainerModel: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.SelectionContainerModel),
/* harmony export */   SelectionModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectionModel),
/* harmony export */   SelectionRangeSliderModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectionRangeSliderModel),
/* harmony export */   SelectionRangeSliderView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectionRangeSliderView),
/* harmony export */   SelectionSliderModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectionSliderModel),
/* harmony export */   SelectionSliderView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectionSliderView),
/* harmony export */   SelectionView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.SelectionView),
/* harmony export */   SliderStyleModel: () => (/* reexport safe */ _widget_int__WEBPACK_IMPORTED_MODULE_13__.SliderStyleModel),
/* harmony export */   StackModel: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.StackModel),
/* harmony export */   StackView: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.StackView),
/* harmony export */   StringModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.StringModel),
/* harmony export */   StringView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.StringView),
/* harmony export */   TabModel: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.TabModel),
/* harmony export */   TabView: () => (/* reexport safe */ _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__.TabView),
/* harmony export */   TagsInputModel: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.TagsInputModel),
/* harmony export */   TagsInputView: () => (/* reexport safe */ _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__.TagsInputView),
/* harmony export */   TextModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.TextModel),
/* harmony export */   TextStyleModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.TextStyleModel),
/* harmony export */   TextView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.TextView),
/* harmony export */   TextareaModel: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.TextareaModel),
/* harmony export */   TextareaView: () => (/* reexport safe */ _widget_string__WEBPACK_IMPORTED_MODULE_19__.TextareaView),
/* harmony export */   TimeModel: () => (/* reexport safe */ _widget_time__WEBPACK_IMPORTED_MODULE_12__.TimeModel),
/* harmony export */   TimeView: () => (/* reexport safe */ _widget_time__WEBPACK_IMPORTED_MODULE_12__.TimeView),
/* harmony export */   ToggleButtonModel: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.ToggleButtonModel),
/* harmony export */   ToggleButtonStyleModel: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.ToggleButtonStyleModel),
/* harmony export */   ToggleButtonView: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.ToggleButtonView),
/* harmony export */   ToggleButtonsModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.ToggleButtonsModel),
/* harmony export */   ToggleButtonsStyleModel: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.ToggleButtonsStyleModel),
/* harmony export */   ToggleButtonsView: () => (/* reexport safe */ _widget_selection__WEBPACK_IMPORTED_MODULE_16__.ToggleButtonsView),
/* harmony export */   VBoxModel: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.VBoxModel),
/* harmony export */   VBoxView: () => (/* reexport safe */ _widget_box__WEBPACK_IMPORTED_MODULE_5__.VBoxView),
/* harmony export */   ValidModel: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.ValidModel),
/* harmony export */   ValidView: () => (/* reexport safe */ _widget_bool__WEBPACK_IMPORTED_MODULE_3__.ValidView),
/* harmony export */   VideoModel: () => (/* reexport safe */ _widget_video__WEBPACK_IMPORTED_MODULE_7__.VideoModel),
/* harmony export */   VideoView: () => (/* reexport safe */ _widget_video__WEBPACK_IMPORTED_MODULE_7__.VideoView),
/* harmony export */   datetime_serializers: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.datetime_serializers),
/* harmony export */   deserialize_date: () => (/* reexport safe */ _widget_date__WEBPACK_IMPORTED_MODULE_10__.deserialize_date),
/* harmony export */   deserialize_datetime: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.deserialize_datetime),
/* harmony export */   deserialize_naive: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.deserialize_naive),
/* harmony export */   deserialize_time: () => (/* reexport safe */ _widget_time__WEBPACK_IMPORTED_MODULE_12__.deserialize_time),
/* harmony export */   escape_html: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_0__.escape_html),
/* harmony export */   naive_serializers: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.naive_serializers),
/* harmony export */   reject: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_0__.reject),
/* harmony export */   resolvePromisesDict: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_0__.resolvePromisesDict),
/* harmony export */   serialize_date: () => (/* reexport safe */ _widget_date__WEBPACK_IMPORTED_MODULE_10__.serialize_date),
/* harmony export */   serialize_datetime: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.serialize_datetime),
/* harmony export */   serialize_naive: () => (/* reexport safe */ _widget_datetime__WEBPACK_IMPORTED_MODULE_11__.serialize_naive),
/* harmony export */   serialize_time: () => (/* reexport safe */ _widget_time__WEBPACK_IMPORTED_MODULE_12__.serialize_time),
/* harmony export */   time_serializers: () => (/* reexport safe */ _widget_time__WEBPACK_IMPORTED_MODULE_12__.time_serializers),
/* harmony export */   typeset: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_0__.typeset),
/* harmony export */   uuid: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_0__.uuid),
/* harmony export */   version: () => (/* binding */ version)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./version */ "../../packages/controls/lib/version.js");
/* harmony import */ var _widget_link__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./widget_link */ "../../packages/controls/lib/widget_link.js");
/* harmony import */ var _widget_bool__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./widget_bool */ "../../packages/controls/lib/widget_bool.js");
/* harmony import */ var _widget_button__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./widget_button */ "../../packages/controls/lib/widget_button.js");
/* harmony import */ var _widget_box__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./widget_box */ "../../packages/controls/lib/widget_box.js");
/* harmony import */ var _widget_image__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./widget_image */ "../../packages/controls/lib/widget_image.js");
/* harmony import */ var _widget_video__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./widget_video */ "../../packages/controls/lib/widget_video.js");
/* harmony import */ var _widget_audio__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./widget_audio */ "../../packages/controls/lib/widget_audio.js");
/* harmony import */ var _widget_color__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./widget_color */ "../../packages/controls/lib/widget_color.js");
/* harmony import */ var _widget_date__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(/*! ./widget_date */ "../../packages/controls/lib/widget_date.js");
/* harmony import */ var _widget_datetime__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(/*! ./widget_datetime */ "../../packages/controls/lib/widget_datetime.js");
/* harmony import */ var _widget_time__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(/*! ./widget_time */ "../../packages/controls/lib/widget_time.js");
/* harmony import */ var _widget_int__WEBPACK_IMPORTED_MODULE_13__ = __webpack_require__(/*! ./widget_int */ "../../packages/controls/lib/widget_int.js");
/* harmony import */ var _widget_float__WEBPACK_IMPORTED_MODULE_14__ = __webpack_require__(/*! ./widget_float */ "../../packages/controls/lib/widget_float.js");
/* harmony import */ var _widget_controller__WEBPACK_IMPORTED_MODULE_15__ = __webpack_require__(/*! ./widget_controller */ "../../packages/controls/lib/widget_controller.js");
/* harmony import */ var _widget_selection__WEBPACK_IMPORTED_MODULE_16__ = __webpack_require__(/*! ./widget_selection */ "../../packages/controls/lib/widget_selection.js");
/* harmony import */ var _widget_selectioncontainer__WEBPACK_IMPORTED_MODULE_17__ = __webpack_require__(/*! ./widget_selectioncontainer */ "../../packages/controls/lib/widget_selectioncontainer.js");
/* harmony import */ var _widget_tagsinput__WEBPACK_IMPORTED_MODULE_18__ = __webpack_require__(/*! ./widget_tagsinput */ "../../packages/controls/lib/widget_tagsinput.js");
/* harmony import */ var _widget_string__WEBPACK_IMPORTED_MODULE_19__ = __webpack_require__(/*! ./widget_string */ "../../packages/controls/lib/widget_string.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_20__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _widget_upload__WEBPACK_IMPORTED_MODULE_21__ = __webpack_require__(/*! ./widget_upload */ "../../packages/controls/lib/widget_upload.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.






















const version = (__webpack_require__(/*! ../package.json */ "../../packages/controls/package.json").version);


/***/ }),

/***/ "../../packages/controls/lib/lumino/accordion.js":
/*!*******************************************************!*\
  !*** ../../packages/controls/lib/lumino/accordion.js ***!
  \*******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Accordion: () => (/* binding */ Accordion),
/* harmony export */   Collapse: () => (/* binding */ Collapse)
/* harmony export */ });
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @lumino/algorithm */ "webpack/sharing/consume/default/@lumino/algorithm");
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @lumino/signaling */ "webpack/sharing/consume/default/@lumino/signaling");
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_lumino_signaling__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @lumino/widgets */ "webpack/sharing/consume/default/@lumino/widgets");
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_lumino_widgets__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _currentselection__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./currentselection */ "../../packages/controls/lib/lumino/currentselection.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.




/**
 * The class name added to Collapse instances.
 */
const COLLAPSE_CLASS = 'jupyter-widget-Collapse';
/**
 * The class name added to a Collapse's header.
 */
const COLLAPSE_HEADER_CLASS = 'jupyter-widget-Collapse-header';
/**
 * The class name added to a Collapse's contents.
 */
const COLLAPSE_CONTENTS_CLASS = 'jupyter-widget-Collapse-contents';
/**
 * The class name added to a Collapse when it is opened
 */
const COLLAPSE_CLASS_OPEN = 'jupyter-widget-Collapse-open';
/**
 * A panel that supports a collapsible header, made from the widget's title.
 * Clicking on the title expands or contracts the widget.
 */
class Collapse extends _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Widget {
    constructor(options) {
        super(options);
        this._collapseChanged = new _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__.Signal(this);
        this.addClass(COLLAPSE_CLASS);
        this._header = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Widget();
        this._header.addClass(COLLAPSE_HEADER_CLASS);
        this._header.node.addEventListener('click', this);
        // Fontawesome icon for caret
        const icon = document.createElement('i');
        icon.classList.add('fa', 'fa-fw', 'fa-caret-right');
        this._header.node.appendChild(icon);
        // Label content
        this._header.node.appendChild(document.createElement('span'));
        this._content = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Panel();
        this._content.addClass(COLLAPSE_CONTENTS_CLASS);
        const layout = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.PanelLayout();
        this.layout = layout;
        layout.addWidget(this._header);
        layout.addWidget(this._content);
        if (options.widget) {
            this.widget = options.widget;
        }
        this.collapsed = false;
    }
    dispose() {
        if (this.isDisposed) {
            return;
        }
        super.dispose();
        this._header = null;
        this._widget = null;
        this._content = null;
    }
    get widget() {
        return this._widget;
    }
    set widget(widget) {
        const oldWidget = this._widget;
        if (oldWidget) {
            oldWidget.disposed.disconnect(this._onChildDisposed, this);
            oldWidget.title.changed.disconnect(this._onTitleChanged, this);
            oldWidget.parent = null;
        }
        this._widget = widget;
        widget.disposed.connect(this._onChildDisposed, this);
        widget.title.changed.connect(this._onTitleChanged, this);
        this._onTitleChanged(widget.title);
        this._content.addWidget(widget);
    }
    get collapsed() {
        return this._collapsed;
    }
    set collapsed(value) {
        // TODO: should we have this check here?
        if (value === this._collapsed) {
            return;
        }
        if (value) {
            this._collapse();
        }
        else {
            this._uncollapse();
        }
    }
    toggle() {
        this.collapsed = !this.collapsed;
    }
    get collapseChanged() {
        return this._collapseChanged;
    }
    _collapse() {
        this._collapsed = true;
        if (this._content) {
            this._content.hide();
        }
        this.removeClass(COLLAPSE_CLASS_OPEN);
        this._header.node.children[0].classList.add('fa-caret-right');
        this._header.node.children[0].classList.remove('fa-caret-down');
        this._collapseChanged.emit(void 0);
    }
    _uncollapse() {
        this._collapsed = false;
        if (this._content) {
            this._content.show();
        }
        this.addClass(COLLAPSE_CLASS_OPEN);
        this._header.node.children[0].classList.add('fa-caret-down');
        this._header.node.children[0].classList.remove('fa-caret-right');
        this._collapseChanged.emit(void 0);
    }
    /**
     * Handle the DOM events for the Collapse widget.
     *
     * @param event - The DOM event sent to the panel.
     *
     * #### Notes
     * This method implements the DOM `EventListener` interface and is
     * called in response to events on the panel's DOM node. It should
     * not be called directly by user code.
     */
    handleEvent(event) {
        switch (event.type) {
            case 'click':
                this._evtClick(event);
                break;
            default:
                break;
        }
    }
    _evtClick(event) {
        this.toggle();
    }
    /**
     * Handle the `changed` signal of a title object.
     */
    _onTitleChanged(sender) {
        this._header.node.children[1].textContent = this._widget.title.label;
    }
    _onChildDisposed(sender) {
        this.dispose();
    }
}
/**
 * The class name added to Accordion instances.
 */
const ACCORDION_CLASS = 'jupyter-widget-Accordion';
/**
 * The class name added to an Accordion child.
 */
const ACCORDION_CHILD_CLASS = 'jupyter-widget-Accordion-child';
const ACCORDION_CHILD_ACTIVE_CLASS = 'jupyter-widget-Accordion-child-active';
/**
 * A panel that supports a collapsible header, made from the widget's title.
 * Clicking on the title expands or contracts the widget.
 */
class Accordion extends _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Panel {
    constructor(options) {
        super(options);
        this._selection = new _currentselection__WEBPACK_IMPORTED_MODULE_3__.Selection(this.widgets);
        this._selection.selectionChanged.connect(this._onSelectionChanged, this);
        this.addClass(ACCORDION_CLASS);
    }
    /**
     * A read-only sequence of the widgets in the panel.
     *
     * #### Notes
     * This is a read-only property.
     */
    /*  get widgets(): ISequence<Widget> {
      return new ArraySequence(toArray(map((this.layout as PanelLayout).widgets, (w: Collapse) => w.widget)));
    }
  */
    get collapseWidgets() {
        return this.layout.widgets;
    }
    get selection() {
        return this._selection;
    }
    indexOf(widget) {
        return _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__.ArrayExt.findFirstIndex(this.collapseWidgets, (w) => w.widget === widget);
    }
    /**
     * Add a widget to the end of the accordion.
     *
     * @param widget - The widget to add to the accordion.
     *
     * @returns The Collapse widget wrapping the added widget.
     *
     * #### Notes
     * The widget will be wrapped in a CollapsedWidget.
     */
    addWidget(widget) {
        const collapse = this._wrapWidget(widget);
        collapse.collapsed = true;
        super.addWidget(collapse);
        this._selection.adjustSelectionForInsert(this.widgets.length - 1, collapse);
        return collapse;
    }
    /**
     * Insert a widget at the specified index.
     *
     * @param index - The index at which to insert the widget.
     *
     * @param widget - The widget to insert into to the accordion.
     *
     * #### Notes
     * If the widget is already contained in the panel, it will be moved.
     */
    insertWidget(index, widget) {
        const collapse = this._wrapWidget(widget);
        collapse.collapsed = true;
        super.insertWidget(index, collapse);
        this._selection.adjustSelectionForInsert(index, collapse);
    }
    removeWidget(widget) {
        const index = this.indexOf(widget);
        if (index >= 0) {
            const collapse = this.collapseWidgets[index];
            widget.parent = null;
            collapse.dispose();
            this._selection.adjustSelectionForRemove(index, null);
        }
    }
    _wrapWidget(widget) {
        const collapse = new Collapse({ widget });
        collapse.addClass(ACCORDION_CHILD_CLASS);
        collapse.collapseChanged.connect(this._onCollapseChange, this);
        return collapse;
    }
    _onCollapseChange(sender) {
        if (!sender.collapsed) {
            this._selection.value = sender;
        }
        else if (this._selection.value === sender && sender.collapsed) {
            this._selection.value = null;
        }
    }
    _onSelectionChanged(sender, change) {
        // Collapse previous widget, open current widget
        const pv = change.previousValue;
        const cv = change.currentValue;
        if (pv) {
            pv.collapsed = true;
            pv.removeClass(ACCORDION_CHILD_ACTIVE_CLASS);
        }
        if (cv) {
            cv.collapsed = false;
            cv.addClass(ACCORDION_CHILD_ACTIVE_CLASS);
        }
    }
}


/***/ }),

/***/ "../../packages/controls/lib/lumino/currentselection.js":
/*!**************************************************************!*\
  !*** ../../packages/controls/lib/lumino/currentselection.js ***!
  \**************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Selection: () => (/* binding */ Selection)
/* harmony export */ });
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @lumino/algorithm */ "webpack/sharing/consume/default/@lumino/algorithm");
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @lumino/signaling */ "webpack/sharing/consume/default/@lumino/signaling");
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_lumino_signaling__WEBPACK_IMPORTED_MODULE_1__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.
/**
 * A variety of convenience methods for maintaining a current selection
 */


class Selection {
    constructor(sequence, options = {}) {
        this._array = null;
        this._value = null;
        this._previousValue = null;
        this._selectionChanged = new _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__.Signal(this);
        this._array = sequence;
        this._insertBehavior = options.insertBehavior || 'select-item-if-needed';
        this._removeBehavior = options.removeBehavior || 'select-item-after';
    }
    /**
     * A signal emitted when the current item is changed.
     *
     * #### Notes
     * This signal is emitted when the currently selected item is changed either
     * through user or programmatic interaction.
     *
     * Notably, this signal is not emitted when the index of the current item
     * changes due to other items being inserted, removed, or moved, but the
     * current item remains the same. It is only emitted when the actual current
     * item is changed.
     */
    get selectionChanged() {
        return this._selectionChanged;
    }
    /**
     * Adjust for setting an item.
     *
     * This should be called *after* the set.
     *
     * @param index - The index set.
     * @param oldValue - The old value at the index.
     */
    adjustSelectionForSet(index) {
        // We just need to send a signal if the currentValue changed.
        // Get the current index and value.
        const pi = this.index;
        const pv = this.value;
        // Exit early if this doesn't affect the selection
        if (index !== pi) {
            return;
        }
        this._updateSelectedValue();
        const cv = this.value;
        // The previous item is now null, since it is no longer in the array.
        this._previousValue = null;
        // Send signal if there was a change
        if (pv !== cv) {
            // Emit the current changed signal.
            this._selectionChanged.emit({
                previousIndex: pi,
                previousValue: pv,
                currentIndex: pi,
                currentValue: cv,
            });
        }
    }
    /**
     * Get the currently selected item.
     *
     * #### Notes
     * This will be `null` if no item is selected.
     */
    get value() {
        return this._value;
    }
    /**
     * Set the currently selected item.
     *
     * #### Notes
     * If the item does not exist in the vector, the currentValue will be set to
     * `null`. This selects the first entry equal to the desired item.
     */
    set value(value) {
        if (value === null || this._array === null) {
            this.index = null;
        }
        else {
            this.index = _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__.ArrayExt.firstIndexOf(this._array, value);
        }
    }
    /**
     * Get the index of the currently selected item.
     *
     * #### Notes
     * This will be `null` if no item is selected.
     */
    get index() {
        return this._index;
    }
    /**
     * Set the index of the currently selected tab.
     *
     * @param index - The index to select.
     *
     * #### Notes
     * If the value is out of range, the index will be set to `null`, which
     * indicates no item is selected.
     */
    set index(index) {
        // Coerce the value to an index.
        let i;
        if (index !== null && this._array !== null) {
            i = Math.floor(index);
            if (i < 0 || i >= this._array.length) {
                i = null;
            }
        }
        else {
            i = null;
        }
        // Bail early if the index will not change.
        if (this._index === i) {
            return;
        }
        // Look up the previous index and item.
        const pi = this._index;
        const pv = this._value;
        // Update the state
        this._index = i;
        this._updateSelectedValue();
        this._previousValue = pv;
        // Emit the current changed signal.
        this._selectionChanged.emit({
            previousIndex: pi,
            previousValue: pv,
            currentIndex: i,
            currentValue: this._value,
        });
    }
    /**
     * Get the selection behavior when inserting a tab.
     */
    get insertBehavior() {
        return this._insertBehavior;
    }
    /**
     * Set the selection behavior when inserting a tab.
     */
    set insertBehavior(value) {
        this._insertBehavior = value;
    }
    /**
     * Get the selection behavior when removing a tab.
     */
    get removeBehavior() {
        return this._removeBehavior;
    }
    /**
     * Set the selection behavior when removing a tab.
     */
    set removeBehavior(value) {
        this._removeBehavior = value;
    }
    /**
     * Adjust the current index for a tab insert operation.
     *
     * @param i - The new index of the inserted item.
     * @param j - The inserted item.
     *
     * #### Notes
     * This method accounts for the tab bar's insertion behavior when adjusting
     * the current index and emitting the changed signal. This should be called
     * after the insertion.
     */
    adjustSelectionForInsert(i, item) {
        // Lookup commonly used variables.
        const cv = this._value;
        const ci = this._index;
        const bh = this._insertBehavior;
        // Handle the behavior where the new item is always selected,
        // or the behavior where the new item is selected if needed.
        if (bh === 'select-item' ||
            (bh === 'select-item-if-needed' && ci === null)) {
            this._index = i;
            this._value = item;
            this._previousValue = cv;
            this._selectionChanged.emit({
                previousIndex: ci,
                previousValue: cv,
                currentIndex: i,
                currentValue: item,
            });
            return;
        }
        // Otherwise, silently adjust the current index if needed.
        if (ci !== null && ci >= i) {
            this._index++;
        }
    }
    /**
     * Clear the selection and history.
     */
    clearSelection() {
        // Get the current index and item.
        const pi = this._index;
        const pv = this._value;
        // Reset the current index and previous item.
        this._index = null;
        this._value = null;
        this._previousValue = null;
        // If no item was selected, there's nothing else to do.
        if (pi === null) {
            return;
        }
        // Emit the current changed signal.
        this._selectionChanged.emit({
            previousIndex: pi,
            previousValue: pv,
            currentIndex: this._index,
            currentValue: this._value,
        });
    }
    /**
     * Adjust the current index for an item remove operation.
     *
     * @param i - The former index of the removed item.
     * @param item - The removed item.
     *
     * #### Notes
     * This method accounts for the remove behavior when adjusting the current
     * index and emitting the changed signal. It should be called after the item
     * is removed.
     */
    adjustSelectionForRemove(i, item) {
        // If we have no selection, there is nothing to do
        if (this._index === null) {
            return;
        }
        // Lookup commonly used variables.
        const ci = this._index;
        const bh = this._removeBehavior;
        // Silently adjust the index if the current item is not removed.
        if (ci !== i) {
            if (ci > i) {
                this._index--;
            }
            return;
        }
        // No item gets selected if the vector is empty.
        if (!this._array || this._array.length === 0) {
            // Reset the current index and previous item.
            this._index = null;
            this._value = null;
            this._previousValue = null;
            this._selectionChanged.emit({
                previousIndex: i,
                previousValue: item,
                currentIndex: this._index,
                currentValue: this._value,
            });
            return;
        }
        // Handle behavior where the next sibling item is selected.
        if (bh === 'select-item-after') {
            this._index = Math.min(i, this._array.length - 1);
            this._updateSelectedValue();
            this._previousValue = null;
            this._selectionChanged.emit({
                previousIndex: i,
                previousValue: item,
                currentIndex: this._index,
                currentValue: this._value,
            });
            return;
        }
        // Handle behavior where the previous sibling item is selected.
        if (bh === 'select-item-before') {
            this._index = Math.max(0, i - 1);
            this._updateSelectedValue();
            this._previousValue = null;
            this._selectionChanged.emit({
                previousIndex: i,
                previousValue: item,
                currentIndex: this._index,
                currentValue: this._value,
            });
            return;
        }
        // Handle behavior where the previous history item is selected.
        if (bh === 'select-previous-item') {
            if (this._previousValue) {
                this.value = this._previousValue;
            }
            else {
                this._index = Math.min(i, this._array.length - 1);
                this._updateSelectedValue();
            }
            this._previousValue = null;
            this._selectionChanged.emit({
                previousIndex: i,
                previousValue: item,
                currentIndex: this._index,
                currentValue: this.value,
            });
            return;
        }
        // Otherwise, no item gets selected.
        this._index = null;
        this._value = null;
        this._previousValue = null;
        this._selectionChanged.emit({
            previousIndex: i,
            previousValue: item,
            currentIndex: this._index,
            currentValue: this._value,
        });
    }
    /**
     * Set the current value based on the current index.
     */
    _updateSelectedValue() {
        const i = this._index;
        this._value = i !== null && this._array ? this._array[i] : null;
    }
}


/***/ }),

/***/ "../../packages/controls/lib/lumino/tabpanel.js":
/*!******************************************************!*\
  !*** ../../packages/controls/lib/lumino/tabpanel.js ***!
  \******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   EventedPanel: () => (/* binding */ EventedPanel),
/* harmony export */   TabPanel: () => (/* binding */ TabPanel)
/* harmony export */ });
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @lumino/messaging */ "webpack/sharing/consume/default/@lumino/messaging");
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_lumino_messaging__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @lumino/signaling */ "webpack/sharing/consume/default/@lumino/signaling");
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_lumino_signaling__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _lumino_domutils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @lumino/domutils */ "webpack/sharing/consume/default/@lumino/domutils");
/* harmony import */ var _lumino_domutils__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_lumino_domutils__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @lumino/widgets */ "webpack/sharing/consume/default/@lumino/widgets");
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(_lumino_widgets__WEBPACK_IMPORTED_MODULE_3__);
/* This file has code derived from Lumino. The license for this Lumino code is:

Copyright (c) 2019 Project Jupyter Contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Copyright (c) 2014-2017, PhosphorJS Contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/




/**
 * A panel where visible widgets are stacked atop one another.
 *
 * #### Notes
 * This class provides a convenience wrapper around a [[PanelLayout]].
 */
class EventedPanel extends _lumino_widgets__WEBPACK_IMPORTED_MODULE_3__.Panel {
    constructor() {
        super(...arguments);
        this._widgetRemoved = new _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__.Signal(this);
    }
    /**
     * A signal emitted when a widget is removed from the panel.
     */
    get widgetRemoved() {
        return this._widgetRemoved;
    }
    /**
     * A message handler invoked on a `'child-removed'` message.
     */
    onChildRemoved(msg) {
        this._widgetRemoved.emit(msg.child);
    }
}
/**
 * A widget which combines a `TabBar` and a `EventedPanel`.
 *
 * #### Notes
 * This is a simple panel which handles the common case of a tab bar
 * placed next to a content area. The selected tab controls the widget
 * which is shown in the content area.
 *
 * For use cases which require more control than is provided by this
 * panel, the `TabBar` widget may be used independently.
 *
 * TODO: Support setting the direction??
 */
class TabPanel extends _lumino_widgets__WEBPACK_IMPORTED_MODULE_3__.Widget {
    /**
     * Construct a new tab panel.
     *
     * @param options - The options for initializing the tab panel.
     */
    constructor(options = {}) {
        super();
        this._currentChanged = new _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__.Signal(this);
        this.addClass('jupyter-widget-TabPanel');
        // Create the tab bar and contents panel.
        this.tabBar = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_3__.TabBar(options);
        this.tabBar.addClass('jupyter-widget-TabPanel-tabBar');
        this.tabContents = new EventedPanel();
        this.tabContents.addClass('jupyter-widget-TabPanel-tabContents');
        // Connect the tab bar signal handlers.
        this.tabBar.tabMoved.connect(this._onTabMoved, this);
        this.tabBar.currentChanged.connect(this._onCurrentChanged, this);
        this.tabBar.tabCloseRequested.connect(this._onTabCloseRequested, this);
        this.tabBar.tabActivateRequested.connect(this._onTabActivateRequested, this);
        // Connect the evented panel signal handlers.
        this.tabContents.widgetRemoved.connect(this._onWidgetRemoved, this);
        // Create the layout.
        const layout = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_3__.PanelLayout();
        // Add the child widgets to the layout.
        layout.addWidget(this.tabBar);
        layout.addWidget(this.tabContents);
        // Install the layout on the tab panel.
        this.layout = layout;
    }
    /**
     * A signal emitted when the current tab is changed.
     *
     * #### Notes
     * This signal is emitted when the currently selected tab is changed
     * either through user or programmatic interaction.
     *
     * Notably, this signal is not emitted when the index of the current
     * tab changes due to tabs being inserted, removed, or moved. It is
     * only emitted when the actual current tab node is changed.
     */
    get currentChanged() {
        return this._currentChanged;
    }
    /**
     * Get the index of the currently selected tab.
     *
     * #### Notes
     * This will be `null` if no tab is selected.
     */
    get currentIndex() {
        const currentIndex = this.tabBar.currentIndex;
        // Lumino tab bars have an index of -1 if no tab is selected
        return currentIndex === -1 ? null : currentIndex;
    }
    /**
     * Set the index of the currently selected tab.
     *
     * #### Notes
     * If the index is out of range, it will be set to `null`.
     */
    set currentIndex(value) {
        this.tabBar.currentIndex = value === null ? -1 : value;
    }
    /**
     * Get the currently selected widget.
     *
     * #### Notes
     * This will be `null` if there is no selected tab.
     */
    get currentWidget() {
        const title = this.tabBar.currentTitle;
        return title ? title.owner : null;
    }
    /**
     * Set the currently selected widget.
     *
     * #### Notes
     * If the widget is not in the panel, it will be set to `null`.
     */
    set currentWidget(value) {
        this.tabBar.currentTitle = value ? value.title : null;
    }
    /**
     * Get the whether the tabs are movable by the user.
     *
     * #### Notes
     * Tabs can always be moved programmatically.
     */
    get tabsMovable() {
        return this.tabBar.tabsMovable;
    }
    /**
     * Set the whether the tabs are movable by the user.
     *
     * #### Notes
     * Tabs can always be moved programmatically.
     */
    set tabsMovable(value) {
        this.tabBar.tabsMovable = value;
    }
    /**
     * A read-only array of the widgets in the panel.
     */
    get widgets() {
        return this.tabContents.widgets;
    }
    /**
     * Add a widget to the end of the tab panel.
     *
     * @param widget - The widget to add to the tab panel.
     *
     * #### Notes
     * If the widget is already contained in the panel, it will be moved.
     *
     * The widget's `title` is used to populate the tab.
     */
    addWidget(widget) {
        this.insertWidget(this.widgets.length, widget);
    }
    /**
     * Insert a widget into the tab panel at a specified index.
     *
     * @param index - The index at which to insert the widget.
     *
     * @param widget - The widget to insert into to the tab panel.
     *
     * #### Notes
     * If the widget is already contained in the panel, it will be moved.
     *
     * The widget's `title` is used to populate the tab.
     */
    insertWidget(index, widget) {
        if (widget !== this.currentWidget) {
            widget.hide();
        }
        this.tabContents.insertWidget(index, widget);
        this.tabBar.insertTab(index, widget.title);
    }
    /**
     * Handle the `currentChanged` signal from the tab bar.
     */
    _onCurrentChanged(sender, args) {
        // Extract the previous and current title from the args.
        const { previousIndex, previousTitle, currentIndex, currentTitle } = args;
        // Extract the widgets from the titles.
        const previousWidget = previousTitle ? previousTitle.owner : null;
        const currentWidget = currentTitle ? currentTitle.owner : null;
        // Hide the previous widget.
        if (previousWidget) {
            previousWidget.hide();
        }
        // Show the current widget.
        if (currentWidget) {
            currentWidget.show();
        }
        // Emit the `currentChanged` signal for the tab panel.
        this._currentChanged.emit({
            previousIndex,
            previousWidget,
            currentIndex,
            currentWidget,
        });
        // Flush the message loop on IE and Edge to prevent flicker.
        if (_lumino_domutils__WEBPACK_IMPORTED_MODULE_2__.Platform.IS_EDGE || _lumino_domutils__WEBPACK_IMPORTED_MODULE_2__.Platform.IS_IE) {
            _lumino_messaging__WEBPACK_IMPORTED_MODULE_0__.MessageLoop.flush();
        }
    }
    /**
     * Handle the `tabActivateRequested` signal from the tab bar.
     */
    _onTabActivateRequested(sender, args) {
        args.title.owner.activate();
    }
    /**
     * Handle the `tabCloseRequested` signal from the tab bar.
     */
    _onTabCloseRequested(sender, args) {
        args.title.owner.close();
    }
    /**
     * Handle the `tabMoved` signal from the tab bar.
     */
    _onTabMoved(sender, args) {
        this.tabContents.insertWidget(args.toIndex, args.title.owner);
    }
    /**
     * Handle the `widgetRemoved` signal from the stacked panel.
     */
    _onWidgetRemoved(sender, widget) {
        this.tabBar.removeTab(widget.title);
    }
}


/***/ }),

/***/ "../../packages/controls/lib/utils.js":
/*!********************************************!*\
  !*** ../../packages/controls/lib/utils.js ***!
  \********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   escape_html: () => (/* binding */ escape_html),
/* harmony export */   reject: () => (/* binding */ reject),
/* harmony export */   resolvePromisesDict: () => (/* reexport safe */ _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.resolvePromisesDict),
/* harmony export */   typeset: () => (/* binding */ typeset),
/* harmony export */   uuid: () => (/* reexport safe */ _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.uuid)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

/**
 * Apply MathJax rendering to an element, and optionally set its text.
 *
 * If MathJax is not available, make no changes.
 *
 * Parameters
 * ----------
 * element: Node
 * text: optional string
 */
function typeset(element, text) {
    if (text !== void 0) {
        element.textContent = text;
    }
    if (window.MathJax !== void 0) {
        MathJax.Hub.Queue(['Typeset', MathJax.Hub, element]);
    }
}
/**
 * escape text to HTML
 */
function escape_html(text) {
    const esc = document.createElement('div');
    esc.textContent = text;
    return esc.innerHTML;
}
/**
 * Creates a wrappable Promise rejection function.
 */
function reject(message, log) {
    return function promiseRejection(error) {
        if (log) {
            console.error(new Error(message));
        }
        throw error;
    };
}


/***/ }),

/***/ "../../packages/controls/lib/widget_audio.js":
/*!***************************************************!*\
  !*** ../../packages/controls/lib/widget_audio.js ***!
  \***************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   AudioModel: () => (/* binding */ AudioModel),
/* harmony export */   AudioView: () => (/* binding */ AudioView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


class AudioModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'AudioModel', _view_name: 'AudioView', format: 'mp3', autoplay: true, loop: true, controls: true, value: new DataView(new ArrayBuffer(0)) });
    }
}
AudioModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel.serializers), { value: {
        serialize: (value) => {
            return new DataView(value.buffer.slice(0));
        },
    } });
class AudioView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    render() {
        /**
         * Called when view is rendered.
         */
        super.render();
        this.luminoWidget.addClass('jupyter-widgets');
        this.update(); // Set defaults.
    }
    update() {
        /**
         * Update the contents of this view
         *
         * Called when the model is changed.  The model may have been
         * changed by another view or by a state update from the back-end.
         */
        let url;
        const format = this.model.get('format');
        const value = this.model.get('value');
        if (format !== 'url') {
            const blob = new Blob([value], {
                type: `audio/${this.model.get('format')}`,
            });
            url = URL.createObjectURL(blob);
        }
        else {
            url = new TextDecoder('utf-8').decode(value.buffer);
        }
        // Clean up the old objectURL
        const oldurl = this.el.src;
        this.el.src = url;
        if (oldurl) {
            URL.revokeObjectURL(oldurl);
        }
        // Audio attributes
        this.el.loop = this.model.get('loop');
        this.el.autoplay = this.model.get('autoplay');
        this.el.controls = this.model.get('controls');
        return super.update();
    }
    remove() {
        if (this.el.src) {
            URL.revokeObjectURL(this.el.src);
        }
        super.remove();
    }
    preinitialize() {
        // Must set this before the initialize method creates the element
        this.tagName = 'audio';
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_bool.js":
/*!**************************************************!*\
  !*** ../../packages/controls/lib/widget_bool.js ***!
  \**************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   BoolModel: () => (/* binding */ BoolModel),
/* harmony export */   CheckboxModel: () => (/* binding */ CheckboxModel),
/* harmony export */   CheckboxStyleModel: () => (/* binding */ CheckboxStyleModel),
/* harmony export */   CheckboxView: () => (/* binding */ CheckboxView),
/* harmony export */   ToggleButtonModel: () => (/* binding */ ToggleButtonModel),
/* harmony export */   ToggleButtonStyleModel: () => (/* binding */ ToggleButtonStyleModel),
/* harmony export */   ToggleButtonView: () => (/* binding */ ToggleButtonView),
/* harmony export */   ValidModel: () => (/* binding */ ValidModel),
/* harmony export */   ValidView: () => (/* binding */ ValidView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.



class CheckboxStyleModel extends _widget_description__WEBPACK_IMPORTED_MODULE_2__.DescriptionStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'CheckboxStyleModel' });
    }
}
CheckboxStyleModel.styleProperties = Object.assign(Object.assign({}, _widget_description__WEBPACK_IMPORTED_MODULE_2__.DescriptionStyleModel.styleProperties), { background: {
        selector: '',
        attribute: 'background',
        default: null,
    } });
class ToggleButtonStyleModel extends _widget_description__WEBPACK_IMPORTED_MODULE_2__.DescriptionStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ToggleButtonStyleModel' });
    }
}
ToggleButtonStyleModel.styleProperties = Object.assign(Object.assign({}, _widget_description__WEBPACK_IMPORTED_MODULE_2__.DescriptionStyleModel.styleProperties), { font_family: {
        selector: '',
        attribute: 'font-family',
        default: '',
    }, font_size: {
        selector: '',
        attribute: 'font-size',
        default: '',
    }, font_style: {
        selector: '',
        attribute: 'font-style',
        default: '',
    }, font_variant: {
        selector: '',
        attribute: 'font-variant',
        default: '',
    }, font_weight: {
        selector: '',
        attribute: 'font-weight',
        default: '',
    }, text_color: {
        selector: '',
        attribute: 'color',
        default: '',
    }, text_decoration: {
        selector: '',
        attribute: 'text-decoration',
        default: '',
    } });
class BoolModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { value: false, disabled: false, _model_name: 'BoolModel' });
    }
}
class CheckboxModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { indent: true, style: null, _view_name: 'CheckboxView', _model_name: 'CheckboxModel' });
    }
}
class CheckboxView extends _widget_description__WEBPACK_IMPORTED_MODULE_2__.DescriptionView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-checkbox');
        // adding a zero-width space to the label to help
        // the browser set the baseline correctly
        this.label.innerHTML = '&#8203;';
        // label containing the checkbox and description span
        this.checkboxLabel = document.createElement('label');
        this.checkboxLabel.classList.add('widget-label-basic');
        this.el.appendChild(this.checkboxLabel);
        // checkbox
        this.checkbox = document.createElement('input');
        this.checkbox.setAttribute('type', 'checkbox');
        this.checkboxLabel.appendChild(this.checkbox);
        // span to the right of the checkbox that will render the description
        this.descriptionSpan = document.createElement('span');
        this.checkboxLabel.appendChild(this.descriptionSpan);
        this.listenTo(this.model, 'change:indent', this.updateIndent);
        this.listenTo(this.model, 'change:tabbable', this.updateTabindex);
        this.update(); // Set defaults.
        this.updateDescription();
        this.updateIndent();
        this.updateTabindex();
        this.updateTooltip();
    }
    /**
     * Overridden from super class
     *
     * Update the description span (rather than the label) since
     * we want the description to the right of the checkbox.
     */
    updateDescription() {
        // can be called before the view is fully initialized
        if (this.checkboxLabel == null) {
            return;
        }
        const description = this.model.get('description');
        if (this.model.get('description_allow_html')) {
            this.descriptionSpan.innerHTML =
                this.model.widget_manager.inline_sanitize(description);
        }
        else {
            this.descriptionSpan.textContent = description;
        }
        this.typeset(this.descriptionSpan);
        this.descriptionSpan.title = description;
        this.checkbox.title = description;
    }
    /**
     * Update the visibility of the label in the super class
     * to provide the optional indent.
     */
    updateIndent() {
        const indent = this.model.get('indent');
        this.label.style.display = indent ? '' : 'none';
    }
    updateTabindex() {
        if (!this.checkbox) {
            return; // we might be constructing the parent
        }
        const tabbable = this.model.get('tabbable');
        if (tabbable === true) {
            this.checkbox.setAttribute('tabIndex', '0');
        }
        else if (tabbable === false) {
            this.checkbox.setAttribute('tabIndex', '-1');
        }
        else if (tabbable === null) {
            this.checkbox.removeAttribute('tabIndex');
        }
    }
    updateTooltip() {
        if (!this.checkbox)
            return; // we might be constructing the parent
        const title = this.model.get('tooltip');
        if (!title) {
            this.checkbox.removeAttribute('title');
        }
        else if (this.model.get('description').length === 0) {
            this.checkbox.setAttribute('title', title);
        }
    }
    events() {
        return {
            'click input[type="checkbox"]': '_handle_click',
        };
    }
    /**
     * Handles when the checkbox is clicked.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    _handle_click() {
        const value = this.model.get('value');
        this.model.set('value', !value, { updated_view: this });
        this.touch();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        this.checkbox.checked = this.model.get('value');
        if (options === undefined || options.updated_view != this) {
            this.checkbox.disabled = this.model.get('disabled');
        }
        return super.update();
    }
    /**
     * Handle message sent to the front end.
     *
     * Used to focus or blur the widget.
     */
    handle_message(content) {
        if (content.do == 'focus') {
            this.checkbox.focus();
        }
        else if (content.do == 'blur') {
            this.checkbox.blur();
        }
    }
}
class ToggleButtonModel extends BoolModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'ToggleButtonView', _model_name: 'ToggleButtonModel', tooltip: '', icon: '', button_style: '', style: null });
    }
}
class ToggleButtonView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('jupyter-button');
        this.el.classList.add('widget-toggle-button');
        this.listenTo(this.model, 'change:button_style', this.update_button_style);
        this.listenTo(this.model, 'change:tabbable', this.updateTabindex);
        this.set_button_style();
        this.update(); // Set defaults.
    }
    update_button_style() {
        this.update_mapped_classes(ToggleButtonView.class_map, 'button_style');
    }
    set_button_style() {
        this.set_mapped_classes(ToggleButtonView.class_map, 'button_style');
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if (this.model.get('value')) {
            this.el.classList.add('mod-active');
        }
        else {
            this.el.classList.remove('mod-active');
        }
        if (options === undefined || options.updated_view !== this) {
            this.el.disabled = this.model.get('disabled');
            this.el.setAttribute('tabbable', this.model.get('tabbable'));
            this.el.setAttribute('title', this.model.get('tooltip'));
            const description = this.model.get('description');
            const icon = this.model.get('icon');
            if (description.trim().length === 0 && icon.trim().length === 0) {
                this.el.innerHTML = '&nbsp;'; // Preserve button height
            }
            else {
                this.el.textContent = '';
                if (icon.trim().length) {
                    const i = document.createElement('i');
                    this.el.appendChild(i);
                    i.classList.add('fa');
                    i.classList.add('fa-' + icon);
                }
                this.el.appendChild(document.createTextNode(description));
            }
        }
        this.updateTabindex();
        return super.update();
    }
    events() {
        return {
            // Dictionary of events and their handlers.
            click: '_handle_click',
        };
    }
    /**
     * Handles and validates user input.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    _handle_click(event) {
        event.preventDefault();
        const value = this.model.get('value');
        this.model.set('value', !value, { updated_view: this });
        this.touch();
    }
    preinitialize() {
        // Must set this before the initialize method creates the element
        this.tagName = 'button';
    }
}
ToggleButtonView.class_map = {
    primary: ['mod-primary'],
    success: ['mod-success'],
    info: ['mod-info'],
    warning: ['mod-warning'],
    danger: ['mod-danger'],
};
class ValidModel extends BoolModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { readout: 'Invalid', _view_name: 'ValidView', _model_name: 'ValidModel' });
    }
}
class ValidView extends _widget_description__WEBPACK_IMPORTED_MODULE_2__.DescriptionView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-valid');
        this.el.classList.add('widget-inline-hbox');
        this.icon = document.createElement('i');
        this.icon.classList.add('fa', 'fa-fw');
        this.el.appendChild(this.icon);
        this.readout = document.createElement('span');
        this.readout.classList.add('widget-valid-readout');
        this.readout.classList.add('widget-readout');
        this.el.appendChild(this.readout);
        this.update();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        this.el.classList.remove('mod-valid');
        this.el.classList.remove('mod-invalid');
        this.icon.classList.remove('fa-check');
        this.icon.classList.remove('fa-times');
        this.readout.textContent = this.model.get('readout');
        if (this.model.get('value')) {
            this.el.classList.add('mod-valid');
            this.icon.classList.add('fa-check');
        }
        else {
            this.el.classList.add('mod-invalid');
            this.icon.classList.add('fa-times');
        }
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_box.js":
/*!*************************************************!*\
  !*** ../../packages/controls/lib/widget_box.js ***!
  \*************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   BoxModel: () => (/* binding */ BoxModel),
/* harmony export */   BoxView: () => (/* binding */ BoxView),
/* harmony export */   GridBoxModel: () => (/* binding */ GridBoxModel),
/* harmony export */   GridBoxView: () => (/* binding */ GridBoxView),
/* harmony export */   HBoxModel: () => (/* binding */ HBoxModel),
/* harmony export */   HBoxView: () => (/* binding */ HBoxView),
/* harmony export */   VBoxModel: () => (/* binding */ VBoxModel),
/* harmony export */   VBoxView: () => (/* binding */ VBoxView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @lumino/algorithm */ "webpack/sharing/consume/default/@lumino/algorithm");
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @lumino/messaging */ "webpack/sharing/consume/default/@lumino/messaging");
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(_lumino_messaging__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @lumino/widgets */ "webpack/sharing/consume/default/@lumino/widgets");
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(_lumino_widgets__WEBPACK_IMPORTED_MODULE_4__);
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! jquery */ "webpack/sharing/consume/default/jquery/jquery?123b");
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_5___default = /*#__PURE__*/__webpack_require__.n(jquery__WEBPACK_IMPORTED_MODULE_5__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.






class BoxModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'BoxView', _model_name: 'BoxModel', children: [], box_style: '' });
    }
}
BoxModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel.serializers), { children: { deserialize: _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.unpack_models } });
class HBoxModel extends BoxModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'HBoxView', _model_name: 'HBoxModel' });
    }
}
class VBoxModel extends BoxModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'VBoxView', _model_name: 'VBoxModel' });
    }
}
class BoxView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    _createElement(tagName) {
        this.luminoWidget = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.JupyterLuminoPanelWidget({ view: this });
        return this.luminoWidget.node;
    }
    _setElement(el) {
        if (this.el || el !== this.luminoWidget.node) {
            // Boxes don't allow setting the element beyond the initial creation.
            throw new Error('Cannot reset the DOM element.');
        }
        this.el = this.luminoWidget.node;
        this.$el = jquery__WEBPACK_IMPORTED_MODULE_5___default()(this.luminoWidget.node);
    }
    initialize(parameters) {
        super.initialize(parameters);
        this.children_views = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.ViewList(this.add_child_model, null, this);
        this.listenTo(this.model, 'change:children', this.update_children);
        this.listenTo(this.model, 'change:box_style', this.update_box_style);
        this.luminoWidget.addClass('jupyter-widgets');
        this.luminoWidget.addClass('widget-container');
        this.luminoWidget.addClass('widget-box');
    }
    render() {
        super.render();
        this.update_children();
        this.set_box_style();
    }
    update_children() {
        var _a;
        (_a = this.children_views) === null || _a === void 0 ? void 0 : _a.update(this.model.get('children')).then((views) => {
            // Notify all children that their sizes may have changed.
            views.forEach((view) => {
                _lumino_messaging__WEBPACK_IMPORTED_MODULE_3__.MessageLoop.postMessage(view.luminoWidget, _lumino_widgets__WEBPACK_IMPORTED_MODULE_4__.Widget.ResizeMessage.UnknownSize);
            });
        });
    }
    update_box_style() {
        this.update_mapped_classes(BoxView.class_map, 'box_style');
    }
    set_box_style() {
        this.set_mapped_classes(BoxView.class_map, 'box_style');
    }
    add_child_model(model) {
        // we insert a dummy element so the order is preserved when we add
        // the rendered content later.
        const dummy = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_4__.Widget();
        this.luminoWidget.addWidget(dummy);
        return this.create_child_view(model)
            .then((view) => {
            // replace the dummy widget with the new one.
            const i = _lumino_algorithm__WEBPACK_IMPORTED_MODULE_2__.ArrayExt.firstIndexOf(this.luminoWidget.widgets, dummy);
            this.luminoWidget.insertWidget(i, view.luminoWidget);
            dummy.dispose();
            return view;
        })
            .catch((0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.reject)('Could not add child view to box', true));
    }
    remove() {
        this.children_views = null;
        super.remove();
    }
}
BoxView.class_map = {
    success: ['alert', 'alert-success'],
    info: ['alert', 'alert-info'],
    warning: ['alert', 'alert-warning'],
    danger: ['alert', 'alert-danger'],
};
class HBoxView extends BoxView {
    /**
     * Public constructor
     */
    initialize(parameters) {
        super.initialize(parameters);
        this.luminoWidget.addClass('widget-hbox');
    }
}
class VBoxView extends BoxView {
    /**
     * Public constructor
     */
    initialize(parameters) {
        super.initialize(parameters);
        this.luminoWidget.addClass('widget-vbox');
    }
}
class GridBoxView extends BoxView {
    /**
     * Public constructor
     */
    initialize(parameters) {
        super.initialize(parameters);
        this.luminoWidget.addClass('widget-gridbox');
        // display needn't be set to flex and grid
        this.luminoWidget.removeClass('widget-box');
    }
}
class GridBoxModel extends BoxModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'GridBoxView', _model_name: 'GridBoxModel' });
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_button.js":
/*!****************************************************!*\
  !*** ../../packages/controls/lib/widget_button.js ***!
  \****************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ButtonModel: () => (/* binding */ ButtonModel),
/* harmony export */   ButtonStyleModel: () => (/* binding */ ButtonStyleModel),
/* harmony export */   ButtonView: () => (/* binding */ ButtonView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./version */ "../../packages/controls/lib/version.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.



class ButtonStyleModel extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.StyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ButtonStyleModel', _model_module: '@jupyter-widgets/controls', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION });
    }
}
ButtonStyleModel.styleProperties = {
    button_color: {
        selector: '',
        attribute: 'background-color',
        default: null,
    },
    font_family: {
        selector: '',
        attribute: 'font-family',
        default: '',
    },
    font_size: {
        selector: '',
        attribute: 'font-size',
        default: '',
    },
    font_style: {
        selector: '',
        attribute: 'font-style',
        default: '',
    },
    font_variant: {
        selector: '',
        attribute: 'font-variant',
        default: '',
    },
    font_weight: {
        selector: '',
        attribute: 'font-weight',
        default: '',
    },
    text_color: {
        selector: '',
        attribute: 'color',
        default: '',
    },
    text_decoration: {
        selector: '',
        attribute: 'text-decoration',
        default: '',
    },
};
class ButtonModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { description: '', tooltip: '', disabled: false, icon: '', button_style: '', _view_name: 'ButtonView', _model_name: 'ButtonModel', style: null });
    }
}
class ButtonView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('jupyter-button');
        this.el.classList.add('widget-button');
        this.listenTo(this.model, 'change:button_style', this.update_button_style);
        this.listenTo(this.model, 'change:tabbable', this.updateTabindex);
        this.set_button_style();
        this.update(); // Set defaults.
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        this.el.disabled = this.model.get('disabled');
        this.updateTabindex();
        const tooltip = this.model.get('tooltip');
        const description = this.model.get('description');
        const icon = this.model.get('icon');
        this.el.setAttribute('title', tooltip !== null && tooltip !== void 0 ? tooltip : description);
        if (description.length || icon.length) {
            this.el.textContent = '';
            if (icon.length) {
                const i = document.createElement('i');
                i.classList.add('fa');
                i.classList.add(...icon
                    .split(/[\s]+/)
                    .filter(Boolean)
                    .map((v) => `fa-${v}`));
                if (description.length === 0) {
                    i.classList.add('center');
                }
                this.el.appendChild(i);
            }
            this.el.appendChild(document.createTextNode(description));
        }
        return super.update();
    }
    update_button_style() {
        this.update_mapped_classes(ButtonView.class_map, 'button_style');
    }
    set_button_style() {
        this.set_mapped_classes(ButtonView.class_map, 'button_style');
    }
    /**
     * Dictionary of events and handlers
     */
    events() {
        // TODO: return typing not needed in Typescript later than 1.8.x
        // See http://stackoverflow.com/questions/22077023/why-cant-i-indirectly-return-an-object-literal-to-satisfy-an-index-signature-re and https://github.com/Microsoft/TypeScript/pull/7029
        return { click: '_handle_click' };
    }
    /**
     * Handles when the button is clicked.
     */
    _handle_click(event) {
        event.preventDefault();
        this.send({ event: 'click' });
    }
    preinitialize() {
        // Must set this before the initialize method creates the element
        this.tagName = 'button';
    }
}
ButtonView.class_map = {
    primary: ['mod-primary'],
    success: ['mod-success'],
    info: ['mod-info'],
    warning: ['mod-warning'],
    danger: ['mod-danger'],
};


/***/ }),

/***/ "../../packages/controls/lib/widget_color.js":
/*!***************************************************!*\
  !*** ../../packages/controls/lib/widget_color.js ***!
  \***************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ColorPickerModel: () => (/* binding */ ColorPickerModel),
/* harmony export */   ColorPickerView: () => (/* binding */ ColorPickerView)
/* harmony export */ });
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.



const named_colors = {
    aliceblue: '#f0f8ff',
    antiquewhite: '#faebd7',
    aqua: '#00ffff',
    aquamarine: '#7fffd4',
    azure: '#f0ffff',
    beige: '#f5f5dc',
    bisque: '#ffe4c4',
    black: '#000000',
    blanchedalmond: '#ffebcd',
    blue: '#0000ff',
    blueviolet: '#8a2be2',
    brown: '#a52a2a',
    burlywood: '#deb887',
    cadetblue: '#5f9ea0',
    chartreuse: '#7fff00',
    chocolate: '#d2691e',
    coral: '#ff7f50',
    cornflowerblue: '#6495ed',
    cornsilk: '#fff8dc',
    crimson: '#dc143c',
    cyan: '#00ffff',
    darkblue: '#00008b',
    darkcyan: '#008b8b',
    darkgoldenrod: '#b8860b',
    darkgray: '#a9a9a9',
    darkgrey: '#a9a9a9',
    darkgreen: '#006400',
    darkkhaki: '#bdb76b',
    darkmagenta: '#8b008b',
    darkolivegreen: '#556b2f',
    darkorange: '#ff8c00',
    darkorchid: '#9932cc',
    darkred: '#8b0000',
    darksalmon: '#e9967a',
    darkseagreen: '#8fbc8f',
    darkslateblue: '#483d8b',
    darkslategray: '#2f4f4f',
    darkslategrey: '#2f4f4f',
    darkturquoise: '#00ced1',
    darkviolet: '#9400d3',
    deeppink: '#ff1493',
    deepskyblue: '#00bfff',
    dimgray: '#696969',
    dimgrey: '#696969',
    dodgerblue: '#1e90ff',
    firebrick: '#b22222',
    floralwhite: '#fffaf0',
    forestgreen: '#228b22',
    fuchsia: '#ff00ff',
    gainsboro: '#dcdcdc',
    ghostwhite: '#f8f8ff',
    gold: '#ffd700',
    goldenrod: '#daa520',
    gray: '#808080',
    grey: '#808080',
    green: '#008000',
    greenyellow: '#adff2f',
    honeydew: '#f0fff0',
    hotpink: '#ff69b4',
    indianred: '#cd5c5c',
    indigo: '#4b0082',
    ivory: '#fffff0',
    khaki: '#f0e68c',
    lavender: '#e6e6fa',
    lavenderblush: '#fff0f5',
    lawngreen: '#7cfc00',
    lemonchiffon: '#fffacd',
    lightblue: '#add8e6',
    lightcoral: '#f08080',
    lightcyan: '#e0ffff',
    lightgoldenrodyellow: '#fafad2',
    lightgreen: '#90ee90',
    lightgray: '#d3d3d3',
    lightgrey: '#d3d3d3',
    lightpink: '#ffb6c1',
    lightsalmon: '#ffa07a',
    lightseagreen: '#20b2aa',
    lightskyblue: '#87cefa',
    lightslategray: '#778899',
    lightslategrey: '#778899',
    lightsteelblue: '#b0c4de',
    lightyellow: '#ffffe0',
    lime: '#00ff00',
    limegreen: '#32cd32',
    linen: '#faf0e6',
    magenta: '#ff00ff',
    maroon: '#800000',
    mediumaquamarine: '#66cdaa',
    mediumblue: '#0000cd',
    mediumorchid: '#ba55d3',
    mediumpurple: '#9370db',
    mediumseagreen: '#3cb371',
    mediumslateblue: '#7b68ee',
    mediumspringgreen: '#00fa9a',
    mediumturquoise: '#48d1cc',
    mediumvioletred: '#c71585',
    midnightblue: '#191970',
    mintcream: '#f5fffa',
    mistyrose: '#ffe4e1',
    moccasin: '#ffe4b5',
    navajowhite: '#ffdead',
    navy: '#000080',
    oldlace: '#fdf5e6',
    olive: '#808000',
    olivedrab: '#6b8e23',
    orange: '#ffa500',
    orangered: '#ff4500',
    orchid: '#da70d6',
    palegoldenrod: '#eee8aa',
    palegreen: '#98fb98',
    paleturquoise: '#afeeee',
    palevioletred: '#db7093',
    papayawhip: '#ffefd5',
    peachpuff: '#ffdab9',
    peru: '#cd853f',
    pink: '#ffc0cb',
    plum: '#dda0dd',
    powderblue: '#b0e0e6',
    purple: '#800080',
    red: '#ff0000',
    rosybrown: '#bc8f8f',
    royalblue: '#4169e1',
    saddlebrown: '#8b4513',
    salmon: '#fa8072',
    sandybrown: '#f4a460',
    seagreen: '#2e8b57',
    seashell: '#fff5ee',
    sienna: '#a0522d',
    silver: '#c0c0c0',
    skyblue: '#87ceeb',
    slateblue: '#6a5acd',
    slategray: '#708090',
    slategrey: '#708090',
    snow: '#fffafa',
    springgreen: '#00ff7f',
    steelblue: '#4682b4',
    tan: '#d2b48c',
    teal: '#008080',
    thistle: '#d8bfd8',
    tomato: '#ff6347',
    turquoise: '#40e0d0',
    violet: '#ee82ee',
    wheat: '#f5deb3',
    white: '#ffffff',
    whitesmoke: '#f5f5f5',
    yellow: '#ffff00',
    yellowgreen: '#9acd32',
};
class ColorPickerModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { value: 'black', concise: false, _model_name: 'ColorPickerModel', _view_name: 'ColorPickerView' });
    }
}
class ColorPickerView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-colorpicker');
        this._color_container = document.createElement('div');
        this._color_container.className =
            'widget-inline-hbox widget-colorpicker-input';
        this.el.appendChild(this._color_container);
        this._textbox = document.createElement('input');
        this._textbox.setAttribute('type', 'text');
        this._textbox.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.uuid)();
        this._color_container.appendChild(this._textbox);
        this._textbox.value = this.model.get('value');
        this._colorpicker = document.createElement('input');
        this._colorpicker.setAttribute('type', 'color');
        this._color_container.appendChild(this._colorpicker);
        this.listenTo(this.model, 'change:value', this._update_value);
        this.listenTo(this.model, 'change:concise', this._update_concise);
        this._update_concise();
        this._update_value();
        this.update();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if (options === undefined || options.updated_view != this) {
            const disabled = this.model.get('disabled');
            this._textbox.disabled = disabled;
            this._colorpicker.disabled = disabled;
        }
        return super.update();
    }
    events() {
        // Typescript doesn't understand that these functions are called, so we
        // specifically use them here so it knows they are being used.
        void this._picker_change;
        void this._text_change;
        return {
            'change [type="color"]': '_picker_change',
            'change [type="text"]': '_text_change',
        };
    }
    _update_value() {
        const value = this.model.get('value');
        this._colorpicker.value = color2hex(value);
        this._textbox.value = value;
    }
    _update_concise() {
        const concise = this.model.get('concise');
        if (concise) {
            this.el.classList.add('concise');
            this._textbox.style.display = 'none';
        }
        else {
            this.el.classList.remove('concise');
            this._textbox.style.display = '';
        }
    }
    _picker_change() {
        this.model.set('value', this._colorpicker.value);
        this.touch();
    }
    _text_change() {
        const value = this._validate_color(this._textbox.value, this.model.get('value'));
        this.model.set('value', value);
        this.touch();
    }
    _validate_color(color, fallback) {
        return color.match(/#[a-fA-F0-9]{3}(?:[a-fA-F0-9]{3})?$/) ||
            named_colors[color.toLowerCase()]
            ? color
            : fallback;
    }
}
/*
 * From a valid html color (named color, 6-digits or 3-digits hex format)
 * return a 6-digits hexadecimal color #rrggbb.
 */
function color2hex(color) {
    return named_colors[color.toLowerCase()] || rgb3_to_rgb6(color);
}
function rgb3_to_rgb6(rgb) {
    if (rgb.length === 7) {
        return rgb;
    }
    else {
        return ('#' +
            rgb.charAt(1) +
            rgb.charAt(1) +
            rgb.charAt(2) +
            rgb.charAt(2) +
            rgb.charAt(3) +
            rgb.charAt(3));
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_controller.js":
/*!********************************************************!*\
  !*** ../../packages/controls/lib/widget_controller.js ***!
  \********************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ControllerAxisModel: () => (/* binding */ ControllerAxisModel),
/* harmony export */   ControllerAxisView: () => (/* binding */ ControllerAxisView),
/* harmony export */   ControllerButtonModel: () => (/* binding */ ControllerButtonModel),
/* harmony export */   ControllerButtonView: () => (/* binding */ ControllerButtonView),
/* harmony export */   ControllerModel: () => (/* binding */ ControllerModel),
/* harmony export */   ControllerView: () => (/* binding */ ControllerView)
/* harmony export */ });
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @lumino/widgets */ "webpack/sharing/consume/default/@lumino/widgets");
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_lumino_widgets__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @lumino/algorithm */ "webpack/sharing/consume/default/@lumino/algorithm");
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! jquery */ "webpack/sharing/consume/default/jquery/jquery?123b");
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_5___default = /*#__PURE__*/__webpack_require__.n(jquery__WEBPACK_IMPORTED_MODULE_5__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.






class ControllerButtonModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ControllerButtonModel', _view_name: 'ControllerButtonView', value: 0.0, pressed: false });
    }
}
/**
 * Very simple view for a gamepad button.
 */
class ControllerButtonView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.DOMWidgetView {
    render() {
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-controller-button');
        this.el.style.width = 'fit-content';
        this.support = document.createElement('div');
        this.support.style.position = 'relative';
        this.support.style.margin = '1px';
        this.support.style.width = '16px';
        this.support.style.height = '16px';
        this.support.style.border = '1px solid black';
        this.support.style.background = 'lightgray';
        this.el.appendChild(this.support);
        this.bar = document.createElement('div');
        this.bar.style.position = 'absolute';
        this.bar.style.width = '100%';
        this.bar.style.bottom = '0px';
        this.bar.style.background = 'gray';
        this.support.appendChild(this.bar);
        this.update();
        this.label = document.createElement('div');
        this.label.textContent = this.model.get('description');
        this.label.style.textAlign = 'center';
        this.el.appendChild(this.label);
    }
    update() {
        this.bar.style.height = 100 * this.model.get('value') + '%';
    }
}
class ControllerAxisModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ControllerAxisModel', _view_name: 'ControllerAxisView', value: 0.0 });
    }
}
/**
 * Very simple view for a gamepad axis.
 */
class ControllerAxisView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.DOMWidgetView {
    render() {
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-controller-axis');
        this.el.style.width = '16px';
        this.el.style.padding = '4px';
        this.support = document.createElement('div');
        this.support.style.position = 'relative';
        this.support.style.margin = '1px';
        this.support.style.width = '4px';
        this.support.style.height = '64px';
        this.support.style.border = '1px solid black';
        this.support.style.background = 'lightgray';
        this.bullet = document.createElement('div');
        this.bullet.style.position = 'absolute';
        this.bullet.style.margin = '-3px';
        this.bullet.style.boxSizing = 'unset';
        this.bullet.style.width = '10px';
        this.bullet.style.height = '10px';
        this.bullet.style.background = 'gray';
        this.label = document.createElement('div');
        this.label.textContent = this.model.get('description');
        this.label.style.textAlign = 'center';
        this.support.appendChild(this.bullet);
        this.el.appendChild(this.support);
        this.el.appendChild(this.label);
        this.update();
    }
    update() {
        this.bullet.style.top = 50 * (this.model.get('value') + 1) + '%';
    }
}
class ControllerModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ControllerModel', _view_name: 'ControllerView', index: 0, name: '', mapping: '', connected: false, timestamp: 0, buttons: [], axes: [] });
    }
    initialize(attributes, options) {
        super.initialize(attributes, options);
        if (navigator.getGamepads === void 0) {
            // Checks if the browser supports the gamepad API
            this.readout = 'This browser does not support gamepads.';
            console.error(this.readout);
        }
        else {
            // Start the wait loop, and listen to updates of the only
            // user-provided attribute, the gamepad index.
            this.readout = 'Connect gamepad and press any button.';
            if (this.get('connected')) {
                // No need to re-create Button and Axis widgets, re-use
                // the models provided by the backend which may already
                // be wired to other things.
                this.update_loop();
            }
            else {
                // Wait for a gamepad to be connected.
                this.wait_loop();
            }
        }
    }
    /**
     * Waits for a gamepad to be connected at the provided index.
     * Once one is connected, it will start the update loop, which
     * populates the update of axes and button values.
     */
    wait_loop() {
        const index = this.get('index');
        const pad = navigator.getGamepads()[index];
        if (pad) {
            this.setup(pad).then((controls) => {
                this.set(controls);
                this.save_changes();
                window.requestAnimationFrame(this.update_loop.bind(this));
            });
        }
        else {
            window.requestAnimationFrame(this.wait_loop.bind(this));
        }
    }
    /**
     * Given a native gamepad object, returns a promise for a dictionary of
     * controls, of the form
     * {
     *     buttons: list of Button models,
     *     axes: list of Axis models,
     * }
     */
    setup(pad) {
        // Set up the main gamepad attributes
        this.set({
            name: pad.id,
            mapping: pad.mapping,
            connected: pad.connected,
            timestamp: pad.timestamp,
        });
        // Create buttons and axes. When done, start the update loop
        return _utils__WEBPACK_IMPORTED_MODULE_4__.resolvePromisesDict({
            buttons: Promise.all(pad.buttons.map((btn, index) => {
                return this._create_button_model(index);
            })),
            axes: Promise.all(pad.axes.map((axis, index) => {
                return this._create_axis_model(index);
            })),
        });
    }
    /**
     * Update axes and buttons values, until the gamepad is disconnected.
     * When the gamepad is disconnected, this.reset_gamepad is called.
     */
    update_loop() {
        const index = this.get('index');
        const id = this.get('name');
        const pad = navigator.getGamepads()[index];
        if (pad && index === pad.index && id === pad.id) {
            this.set({
                timestamp: pad.timestamp,
                connected: pad.connected,
            });
            this.save_changes();
            this.get('buttons').forEach(function (model, index) {
                model.set({
                    value: pad.buttons[index].value,
                    pressed: pad.buttons[index].pressed,
                });
                model.save_changes();
            });
            this.get('axes').forEach(function (model, index) {
                model.set('value', pad.axes[index]);
                model.save_changes();
            });
            window.requestAnimationFrame(this.update_loop.bind(this));
        }
        else {
            this.reset_gamepad();
        }
    }
    /**
     * Resets the gamepad attributes, and start the wait_loop.
     */
    reset_gamepad() {
        this.get('buttons').forEach(function (button) {
            button.close();
        });
        this.get('axes').forEach(function (axis) {
            axis.close();
        });
        this.set({
            name: '',
            mapping: '',
            connected: false,
            timestamp: 0.0,
            buttons: [],
            axes: [],
        });
        this.save_changes();
        window.requestAnimationFrame(this.wait_loop.bind(this));
    }
    /**
     * Creates a gamepad button widget.
     */
    _create_button_model(index) {
        return this.widget_manager
            .new_widget({
            model_name: 'ControllerButtonModel',
            model_module: '@jupyter-widgets/controls',
            model_module_version: this.get('_model_module_version'),
            view_name: 'ControllerButtonView',
            view_module: '@jupyter-widgets/controls',
            view_module_version: this.get('_view_module_version'),
        })
            .then(function (model) {
            model.set('description', index);
            return model;
        });
    }
    /**
     * Creates a gamepad axis widget.
     */
    _create_axis_model(index) {
        return this.widget_manager
            .new_widget({
            model_name: 'ControllerAxisModel',
            model_module: '@jupyter-widgets/controls',
            model_module_version: this.get('_model_module_version'),
            view_name: 'ControllerAxisView',
            view_module: '@jupyter-widgets/controls',
            view_module_version: this.get('_view_module_version'),
        })
            .then(function (model) {
            model.set('description', index);
            return model;
        });
    }
}
ControllerModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDOMWidgetModel.serializers), { buttons: { deserialize: _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.unpack_models }, axes: { deserialize: _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.unpack_models } });
/**
 * A simple view for a gamepad.
 */
class ControllerView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.DOMWidgetView {
    _createElement(tagName) {
        this.luminoWidget = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.JupyterLuminoPanelWidget({ view: this });
        return this.luminoWidget.node;
    }
    _setElement(el) {
        if (this.el || el !== this.luminoWidget.node) {
            // Boxes don't allow setting the element beyond the initial creation.
            throw new Error('Cannot reset the DOM element.');
        }
        this.el = this.luminoWidget.node;
        this.$el = jquery__WEBPACK_IMPORTED_MODULE_5___default()(this.luminoWidget.node);
    }
    initialize(parameters) {
        super.initialize(parameters);
        this.button_views = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.ViewList(this.add_button, null, this);
        this.listenTo(this.model, 'change:buttons', (model, value) => {
            this.button_views.update(value);
        });
        this.axis_views = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.ViewList(this.add_axis, null, this);
        this.listenTo(this.model, 'change:axes', (model, value) => {
            this.axis_views.update(value);
        });
        this.listenTo(this.model, 'change:name', this.update_label);
    }
    render() {
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-controller');
        this.label = document.createElement('div');
        this.el.appendChild(this.label);
        this.axis_box = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Panel();
        this.axis_box.node.style.display = 'flex';
        this.luminoWidget.addWidget(this.axis_box);
        this.button_box = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Panel();
        this.button_box.node.style.display = 'flex';
        this.luminoWidget.addWidget(this.button_box);
        this.button_views.update(this.model.get('buttons'));
        this.axis_views.update(this.model.get('axes'));
        this.update_label();
    }
    update_label() {
        this.label.textContent = this.model.get('name') || this.model.readout;
    }
    add_button(model) {
        // we insert a dummy element so the order is preserved when we add
        // the rendered content later.
        const dummy = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Widget();
        this.button_box.addWidget(dummy);
        return this.create_child_view(model)
            .then((view) => {
            // replace the dummy widget with the new one.
            const i = _lumino_algorithm__WEBPACK_IMPORTED_MODULE_3__.ArrayExt.firstIndexOf(this.button_box.widgets, dummy);
            this.button_box.insertWidget(i, view.luminoWidget);
            dummy.dispose();
            return view;
        })
            .catch((0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.reject)('Could not add child button view to controller', true));
    }
    add_axis(model) {
        // we insert a dummy element so the order is preserved when we add
        // the rendered content later.
        const dummy = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_2__.Widget();
        this.axis_box.addWidget(dummy);
        return this.create_child_view(model)
            .then((view) => {
            // replace the dummy widget with the new one.
            const i = _lumino_algorithm__WEBPACK_IMPORTED_MODULE_3__.ArrayExt.firstIndexOf(this.axis_box.widgets, dummy);
            this.axis_box.insertWidget(i, view.luminoWidget);
            dummy.dispose();
            return view;
        })
            .catch((0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.reject)('Could not add child axis view to controller', true));
    }
    remove() {
        super.remove();
        this.button_views.remove();
        this.axis_views.remove();
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_core.js":
/*!**************************************************!*\
  !*** ../../packages/controls/lib/widget_core.js ***!
  \**************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   CoreDOMWidgetModel: () => (/* binding */ CoreDOMWidgetModel),
/* harmony export */   CoreDescriptionModel: () => (/* binding */ CoreDescriptionModel),
/* harmony export */   CoreWidgetModel: () => (/* binding */ CoreWidgetModel)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./version */ "../../packages/controls/lib/version.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.
// widget_core implements some common patterns for the core widget collection
// that are not to be used directly by third-party widget authors.



class CoreWidgetModel extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.WidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'CoreWidgetModel', _view_module: '@jupyter-widgets/controls', _model_module: '@jupyter-widgets/controls', _view_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION, _model_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION });
    }
}
class CoreDOMWidgetModel extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'CoreDOMWidgetModel', _view_module: '@jupyter-widgets/controls', _model_module: '@jupyter-widgets/controls', _view_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION, _model_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION });
    }
}
class CoreDescriptionModel extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'CoreDescriptionModel', _view_module: '@jupyter-widgets/controls', _model_module: '@jupyter-widgets/controls', _view_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION, _model_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION });
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_date.js":
/*!**************************************************!*\
  !*** ../../packages/controls/lib/widget_date.js ***!
  \**************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DatePickerModel: () => (/* binding */ DatePickerModel),
/* harmony export */   DatePickerView: () => (/* binding */ DatePickerView),
/* harmony export */   deserialize_date: () => (/* binding */ deserialize_date),
/* harmony export */   serialize_date: () => (/* binding */ serialize_date)
/* harmony export */ });
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.



function serialize_date(value) {
    if (value === null) {
        return null;
    }
    else {
        return {
            year: value.getUTCFullYear(),
            month: value.getUTCMonth(),
            date: value.getUTCDate(),
        };
    }
}
function deserialize_date(value) {
    if (value === null) {
        return null;
    }
    else {
        const date = new Date();
        date.setUTCFullYear(value.year, value.month, value.date);
        date.setUTCHours(0, 0, 0, 0);
        return date;
    }
}
class DatePickerModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { value: null, _model_name: 'DatePickerModel', _view_name: 'DatePickerView' });
    }
}
DatePickerModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDescriptionModel.serializers), { value: {
        serialize: serialize_date,
        deserialize: deserialize_date,
    } });
class DatePickerView extends _widget_description__WEBPACK_IMPORTED_MODULE_0__.DescriptionView {
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-datepicker');
        this._datepicker = document.createElement('input');
        this._datepicker.setAttribute('type', 'date');
        this._datepicker.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.uuid)();
        this.el.appendChild(this._datepicker);
        this.listenTo(this.model, 'change:value', this._update_value);
        this._update_value();
        this.update();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if (options === undefined || options.updated_view !== this) {
            this._datepicker.disabled = this.model.get('disabled');
        }
        return super.update();
    }
    events() {
        // Typescript doesn't understand that these functions are called, so we
        // specifically use them here so it knows they are being used.
        void this._picker_change;
        void this._picker_focusout;
        return {
            'change [type="date"]': '_picker_change',
            'focusout [type="date"]': '_picker_focusout',
        };
    }
    _update_value() {
        const value = this.model.get('value');
        this._datepicker.valueAsDate = value;
    }
    _picker_change() {
        if (!this._datepicker.validity.badInput) {
            this.model.set('value', this._datepicker.valueAsDate);
            this.touch();
        }
    }
    _picker_focusout() {
        if (this._datepicker.validity.badInput) {
            this.model.set('value', null);
            this.touch();
        }
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_datetime.js":
/*!******************************************************!*\
  !*** ../../packages/controls/lib/widget_datetime.js ***!
  \******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DatetimeModel: () => (/* binding */ DatetimeModel),
/* harmony export */   DatetimeView: () => (/* binding */ DatetimeView),
/* harmony export */   NaiveDatetimeModel: () => (/* binding */ NaiveDatetimeModel),
/* harmony export */   datetime_serializers: () => (/* binding */ datetime_serializers),
/* harmony export */   deserialize_datetime: () => (/* binding */ deserialize_datetime),
/* harmony export */   deserialize_naive: () => (/* binding */ deserialize_naive),
/* harmony export */   naive_serializers: () => (/* binding */ naive_serializers),
/* harmony export */   serialize_datetime: () => (/* binding */ serialize_datetime),
/* harmony export */   serialize_naive: () => (/* binding */ serialize_naive)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _widget_time__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./widget_time */ "../../packages/controls/lib/widget_time.js");
// Copyright (c) Vidar Tonaas Fauske
// Distributed under the terms of the Modified BSD License.




function serialize_datetime(value) {
    if (value === null) {
        return null;
    }
    else {
        return {
            year: value.getUTCFullYear(),
            month: value.getUTCMonth(),
            date: value.getUTCDate(),
            hours: value.getUTCHours(),
            minutes: value.getUTCMinutes(),
            seconds: value.getUTCSeconds(),
            milliseconds: value.getUTCMilliseconds(),
        };
    }
}
function deserialize_datetime(value) {
    if (value === null) {
        return null;
    }
    else {
        const date = new Date();
        date.setUTCFullYear(value.year, value.month, value.date);
        date.setUTCHours(value.hours, value.minutes, value.seconds, value.milliseconds);
        return date;
    }
}
const datetime_serializers = {
    serialize: serialize_datetime,
    deserialize: deserialize_datetime,
};
class DatetimeModel extends _widget_core__WEBPACK_IMPORTED_MODULE_2__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'DatetimeModel', _view_name: 'DatetimeView', value: null, disabled: false, min: null, max: null });
    }
}
DatetimeModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_2__.CoreDescriptionModel.serializers), { value: datetime_serializers, min: datetime_serializers, max: datetime_serializers });
class DatetimeView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-datetimepicker');
        const test = document.createElement('input');
        test.type = 'datetime-local';
        if (test.type === 'text') {
            // No native support, split into date and time input:
            this._datepicker = document.createElement('input');
            this._datepicker.setAttribute('type', 'date');
            this._datepicker.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.uuid)();
            this._timepicker = document.createElement('input');
            this._timepicker.setAttribute('type', 'time');
            this._timepicker.id = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.uuid)();
            this.el.appendChild(this._datepicker);
            this.el.appendChild(this._timepicker);
        }
        else {
            this._datetimepicker = test;
            this._datetimepicker.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.uuid)();
            this.el.appendChild(this._datetimepicker);
        }
        this.listenTo(this.model, 'change:value', this._update_value);
        this.listenTo(this.model, 'change', this.update2);
        this._update_value();
        this.update2();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update2(model, options) {
        if (options === undefined || options.updated_view !== this) {
            const min = this.model.get('min');
            const max = this.model.get('max');
            if (this._datetimepicker) {
                this._datetimepicker.disabled = this.model.get('disabled');
                this._datetimepicker.min = Private.dt_as_dt_string(min);
                this._datetimepicker.max = Private.dt_as_dt_string(max);
            }
            else {
                this._datepicker.disabled = this.model.get('disabled');
                this._datepicker.min = Private.dt_as_date_string(min);
                this._datepicker.max = Private.dt_as_date_string(max);
                this._timepicker.disabled = this.model.get('disabled');
                // Don't support min/max time here.
                // It could be added by enabling min time if value is min date,
                // and enabling max time if value is max date, but leave as
                // TODO for now.
            }
        }
    }
    events() {
        // Typescript doesn't understand that these functions are called, so we
        // specifically use them here so it knows they are being used.
        void this._picker_change;
        void this._picker_focusout;
        return {
            'change [type="date"]': '_picker_change',
            'change [type="time"]': '_picker_change',
            'change [type="datetime-local"]': '_picker_change',
            'focusout [type="date"]': '_picker_focusout',
            'focusout [type="datetime-local"]': '_picker_focusout',
            'focusout [type="time"]': '_picker_focusout',
        };
    }
    _update_value(model, newValue, options) {
        if (options === undefined || options.updated_view !== this) {
            const value = this.model.get('value');
            if (this._datetimepicker) {
                this._datetimepicker.value = Private.dt_as_dt_string(value);
            }
            else {
                this._datepicker.valueAsDate = value;
                this._timepicker.value = Private.dt_as_time_string(value);
            }
        }
    }
    _picker_change() {
        if (this._datetimepicker) {
            if (!this._datetimepicker.validity.badInput) {
                const v = this._datetimepicker.value;
                let date = v ? new Date(v) : null;
                if (date && isNaN(date.valueOf())) {
                    date = null;
                }
                this.model.set('value', date, { updated_view: this });
                this.touch();
            }
        }
        else {
            if (!this._datepicker.validity.badInput &&
                !this._timepicker.validity.badInput) {
                const date = this._datepicker.valueAsDate;
                const time = (0,_widget_time__WEBPACK_IMPORTED_MODULE_3__.serialize_time)(this._timepicker.value);
                if (date !== null && time !== null) {
                    // * Use local time *
                    date.setHours(time.hours, time.minutes, time.seconds, time.milliseconds);
                }
                this.model.set('value', time !== null && date, { updated_view: this });
                this.touch();
            }
        }
    }
    _picker_focusout() {
        const pickers = [this._datetimepicker, this._datepicker, this._timepicker];
        if (pickers.some((p) => p && p.validity.badInput)) {
            this.model.set('value', null);
            this.touch();
        }
    }
}
var Private;
(function (Private) {
    // eslint-disable-next-line no-inner-declarations
    function dt_as_dt_string(value) {
        if (value === null) {
            return '';
        }
        // Replicate `toISOString()` but in local time zone:
        const parts = [];
        parts.push(`${value.getFullYear().toString().padStart(4, '0')}`);
        parts.push(`-${(value.getMonth() + 1).toString().padStart(2, '0')}`);
        parts.push(`-${value.getDate().toString().padStart(2, '0')}`);
        parts.push(`T${value.getHours().toString().padStart(2, '0')}`);
        parts.push(`:${value.getMinutes().toString().padStart(2, '0')}`);
        if (value.getSeconds() > 0 || value.getMilliseconds() > 0) {
            parts.push(`:${value.getSeconds().toString().padStart(2, '0')}`);
            if (value.getMilliseconds() > 0) {
                parts.push(`.${value.getMilliseconds().toString().padStart(3, '0')}`);
            }
        }
        return parts.join('');
    }
    Private.dt_as_dt_string = dt_as_dt_string;
    // eslint-disable-next-line no-inner-declarations
    function dt_as_date_string(value) {
        return value ? dt_as_dt_string(value).split('T', 2)[0] : '';
    }
    Private.dt_as_date_string = dt_as_date_string;
    // eslint-disable-next-line no-inner-declarations
    function dt_as_time_string(value) {
        return value ? dt_as_dt_string(value).split('T', 2)[1] : '';
    }
    Private.dt_as_time_string = dt_as_time_string;
})(Private || (Private = {}));
function serialize_naive(value) {
    if (value === null) {
        return null;
    }
    else {
        return {
            year: value.getFullYear(),
            month: value.getMonth(),
            date: value.getDate(),
            hours: value.getHours(),
            minutes: value.getMinutes(),
            seconds: value.getSeconds(),
            milliseconds: value.getMilliseconds(),
        };
    }
}
function deserialize_naive(value) {
    if (value === null) {
        return null;
    }
    else {
        const date = new Date();
        date.setFullYear(value.year, value.month, value.date);
        date.setHours(value.hours, value.minutes, value.seconds, value.milliseconds);
        return date;
    }
}
const naive_serializers = {
    serialize: serialize_naive,
    deserialize: deserialize_naive,
};
class NaiveDatetimeModel extends DatetimeModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'NaiveDatetimeModel' });
    }
}
NaiveDatetimeModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_2__.CoreDescriptionModel.serializers), { value: naive_serializers, min: naive_serializers, max: naive_serializers });


/***/ }),

/***/ "../../packages/controls/lib/widget_description.js":
/*!*********************************************************!*\
  !*** ../../packages/controls/lib/widget_description.js ***!
  \*********************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DescriptionModel: () => (/* binding */ DescriptionModel),
/* harmony export */   DescriptionStyleModel: () => (/* binding */ DescriptionStyleModel),
/* harmony export */   DescriptionView: () => (/* binding */ DescriptionView),
/* harmony export */   LabeledDOMWidgetModel: () => (/* binding */ LabeledDOMWidgetModel),
/* harmony export */   LabeledDOMWidgetView: () => (/* binding */ LabeledDOMWidgetView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./version */ "../../packages/controls/lib/version.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.



class DescriptionStyleModel extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.StyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'DescriptionStyleModel', _model_module: '@jupyter-widgets/controls', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION });
    }
}
DescriptionStyleModel.styleProperties = {
    description_width: {
        selector: '.widget-label',
        attribute: 'width',
        default: null,
    },
};
class DescriptionModel extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'DescriptionModel', _view_name: 'DescriptionView', _view_module: '@jupyter-widgets/controls', _model_module: '@jupyter-widgets/controls', _view_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION, _model_module_version: _version__WEBPACK_IMPORTED_MODULE_2__.JUPYTER_CONTROLS_VERSION, description: '', description_allow_html: false });
    }
}
class DescriptionView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    render() {
        this.label = document.createElement('label');
        this.el.appendChild(this.label);
        this.label.className = 'widget-label';
        this.label.style.display = 'none';
        this.listenTo(this.model, 'change:description', this.updateDescription);
        this.listenTo(this.model, 'change:description_allow_html', this.updateDescription);
        this.listenTo(this.model, 'change:tabbable', this.updateTabindex);
        this.updateDescription();
        this.updateTabindex();
        this.updateTooltip();
    }
    typeset(element, text) {
        this.displayed.then(() => {
            var _a;
            const widget_manager = this.model.widget_manager;
            const latexTypesetter = (_a = widget_manager._rendermime) === null || _a === void 0 ? void 0 : _a.latexTypesetter;
            if (latexTypesetter) {
                if (text !== void 0) {
                    element.textContent = text;
                }
                latexTypesetter.typeset(element);
            }
            else {
                return (0,_utils__WEBPACK_IMPORTED_MODULE_1__.typeset)(element, text);
            }
        });
    }
    updateDescription() {
        const description = this.model.get('description');
        if (description.length === 0) {
            this.label.style.display = 'none';
        }
        else {
            if (this.model.get('description_allow_html')) {
                this.label.innerHTML =
                    this.model.widget_manager.inline_sanitize(description);
            }
            else {
                this.label.textContent = description;
            }
            this.typeset(this.label);
            this.label.style.display = '';
        }
    }
    updateTooltip() {
        if (!this.label)
            return;
        this.label.title = this.model.get('tooltip');
    }
}
/**
 * For backwards compatibility with jupyter-js-widgets 2.x.
 *
 * Use DescriptionModel instead.
 */
class LabeledDOMWidgetModel extends DescriptionModel {
}
/**
 * For backwards compatibility with jupyter-js-widgets 2.x.
 *
 * Use DescriptionView instead.
 */
class LabeledDOMWidgetView extends DescriptionView {
}


/***/ }),

/***/ "../../packages/controls/lib/widget_float.js":
/*!***************************************************!*\
  !*** ../../packages/controls/lib/widget_float.js ***!
  \***************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   BoundedFloatModel: () => (/* binding */ BoundedFloatModel),
/* harmony export */   BoundedFloatTextModel: () => (/* binding */ BoundedFloatTextModel),
/* harmony export */   FloatLogSliderModel: () => (/* binding */ FloatLogSliderModel),
/* harmony export */   FloatLogSliderView: () => (/* binding */ FloatLogSliderView),
/* harmony export */   FloatModel: () => (/* binding */ FloatModel),
/* harmony export */   FloatProgressModel: () => (/* binding */ FloatProgressModel),
/* harmony export */   FloatRangeSliderModel: () => (/* binding */ FloatRangeSliderModel),
/* harmony export */   FloatRangeSliderView: () => (/* binding */ FloatRangeSliderView),
/* harmony export */   FloatSliderModel: () => (/* binding */ FloatSliderModel),
/* harmony export */   FloatSliderView: () => (/* binding */ FloatSliderView),
/* harmony export */   FloatTextModel: () => (/* binding */ FloatTextModel),
/* harmony export */   FloatTextView: () => (/* binding */ FloatTextView)
/* harmony export */ });
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _widget_int__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_int */ "../../packages/controls/lib/widget_int.js");
/* harmony import */ var d3_format__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! d3-format */ "../../node_modules/d3-format/src/defaultLocale.js");
/* harmony import */ var nouislider__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! nouislider */ "../../node_modules/nouislider/dist/nouislider.js");
/* harmony import */ var nouislider__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(nouislider__WEBPACK_IMPORTED_MODULE_2__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.




class FloatModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'FloatModel', value: 0 });
    }
}
class BoundedFloatModel extends FloatModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'BoundedFloatModel', max: 100.0, min: 0.0 });
    }
}
class FloatSliderModel extends BoundedFloatModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'FloatSliderModel', _view_name: 'FloatSliderView', step: 1.0, orientation: 'horizontal', _range: false, readout: true, readout_format: '.2f', slider_color: null, continuous_update: true, disabled: false });
    }
    initialize(attributes, options) {
        super.initialize(attributes, options);
        this.on('change:readout_format', this.update_readout_format, this);
        this.update_readout_format();
    }
    update_readout_format() {
        this.readout_formatter = (0,d3_format__WEBPACK_IMPORTED_MODULE_3__.format)(this.get('readout_format'));
    }
}
class FloatLogSliderModel extends BoundedFloatModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'FloatLogSliderModel', _view_name: 'FloatLogSliderView', step: 0.1, orientation: 'horizontal', _range: false, readout: true, readout_format: '.3g', slider_color: null, continuous_update: true, disabled: false, base: 10, value: 1.0, min: 0, max: 4 });
    }
    initialize(attributes, options) {
        super.initialize(attributes, options);
        this.on('change:readout_format', this.update_readout_format, this);
        this.update_readout_format();
    }
    update_readout_format() {
        this.readout_formatter = (0,d3_format__WEBPACK_IMPORTED_MODULE_3__.format)(this.get('readout_format'));
    }
}
class FloatRangeSliderModel extends FloatSliderModel {
}
class FloatSliderView extends _widget_int__WEBPACK_IMPORTED_MODULE_1__.IntSliderView {
    constructor() {
        super(...arguments);
        this._parse_value = parseFloat;
    }
    /**
     * Validate the value of the slider before sending it to the back-end
     * and applying it to the other views on the page.
     */
    _validate_slide_value(x) {
        return x;
    }
}
class FloatLogSliderView extends _widget_int__WEBPACK_IMPORTED_MODULE_1__.BaseIntSliderView {
    constructor() {
        super(...arguments);
        this._parse_value = parseFloat;
    }
    update(options) {
        super.update(options);
        const value = this.model.get('value');
        this.readout.textContent = this.valueToString(value);
    }
    /**
     * Convert from value to exponent
     *
     * @param value the widget value
     * @returns the log-value between the min/max exponents
     */
    logCalc(value) {
        const min = this.model.get('min');
        const max = this.model.get('max');
        const base = this.model.get('base');
        let log_value = Math.log(value) / Math.log(base);
        if (log_value > max) {
            log_value = max;
        }
        else if (log_value < min) {
            log_value = min;
        }
        return log_value;
    }
    createSlider() {
        var _a;
        const orientation = this.model.get('orientation');
        const behavior = this.model.get('behavior');
        nouislider__WEBPACK_IMPORTED_MODULE_2___default().create(this.$slider, {
            start: this.logCalc(this.model.get('value')),
            behaviour: behavior,
            range: {
                min: this.model.get('min'),
                max: this.model.get('max'),
            },
            step: (_a = this.model.get('step')) !== null && _a !== void 0 ? _a : undefined,
            animate: false,
            orientation: orientation,
            direction: orientation === 'horizontal' ? 'ltr' : 'rtl',
            format: {
                from: (value) => Number(value),
                to: (value) => value,
            },
        });
        // Using noUiSlider's 'update' and 'change' events.
        // See reference: https://refreshless.com/nouislider/events-callbacks/
        this.$slider.noUiSlider.on('update', (values, handle) => {
            this.handleSliderUpdateEvent(values, handle);
        });
        this.$slider.noUiSlider.on('change', (values, handle) => {
            this.handleSliderChangeEvent(values, handle);
        });
    }
    /**
     * Write value to a string
     */
    valueToString(value) {
        const format = this.model.readout_formatter;
        return format(value);
    }
    /**
     * Parse value from a string
     */
    stringToValue(text) {
        return text === null ? NaN : this._parse_value(text);
    }
    /**
     * this handles the entry of text into the contentEditable label first, the
     * value is checked if it contains a parseable value then it is clamped
     * within the min-max range of the slider finally, the model is updated if
     * the value is to be changed
     *
     * if any of these conditions are not met, the text is reset
     */
    handleTextChange() {
        let value = this.stringToValue(this.readout.textContent);
        const vmin = this.model.get('min');
        const vmax = this.model.get('max');
        const base = this.model.get('base');
        if (isNaN(value)) {
            this.readout.textContent = this.valueToString(this.model.get('value'));
        }
        else {
            value = Math.max(Math.min(value, Math.pow(base, vmax)), Math.pow(base, vmin));
            if (value !== this.model.get('value')) {
                this.readout.textContent = this.valueToString(value);
                this.model.set('value', value);
                this.touch();
            }
            else {
                this.readout.textContent = this.valueToString(this.model.get('value'));
            }
        }
    }
    /**
     * Called whilst the slider is dragged, tapped or moved by the arrow keys.
     */
    handleSliderUpdateEvent(values, handle) {
        const base = this.model.get('base');
        const actual_value = Math.pow(base, this._validate_slide_value(values[0]));
        this.readout.textContent = this.valueToString(actual_value);
        // Only persist the value while sliding if the continuous_update
        // trait is set to true.
        if (this.model.get('continuous_update')) {
            this.handleSliderChanged(values, handle);
        }
    }
    /**
     * Called when the slider handle is released after dragging,
     * or by tapping or moving by the arrow keys.
     */
    handleSliderChangeEvent(values, handle) {
        const base = this.model.get('base');
        const actual_value = Math.pow(base, this._validate_slide_value(values[0]));
        this.readout.textContent = this.valueToString(actual_value);
        this.handleSliderChanged(values, handle);
    }
    /**
     * Called when the slider value has changed.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    handleSliderChanged(values, handle) {
        if (this._updating_slider) {
            return;
        }
        const base = this.model.get('base');
        const actual_value = Math.pow(base, this._validate_slide_value(values[0]));
        this.model.set('value', actual_value, { updated_view: this });
        this.touch();
    }
    updateSliderValue(model, value, options) {
        if (options.updated_view === this) {
            return;
        }
        const log_value = this.logCalc(this.model.get('value'));
        this.$slider.noUiSlider.set(log_value);
    }
    updateSliderOptions(e) {
        this.$slider.noUiSlider.updateOptions({
            start: this.logCalc(this.model.get('value')),
            range: {
                min: this.model.get('min'),
                max: this.model.get('max'),
            },
            step: this.model.get('step'),
        });
    }
    _validate_slide_value(x) {
        return x;
    }
}
class FloatRangeSliderView extends _widget_int__WEBPACK_IMPORTED_MODULE_1__.IntRangeSliderView {
    constructor() {
        super(...arguments);
        this._parse_value = parseFloat;
        // matches: whitespace?, float, whitespace?, (hyphen, colon, or en-dash), whitespace?, float
        this._range_regex = /^\s*([+-]?(?:\d*\.?\d+|\d+\.)(?:[eE][-:]?\d+)?)\s*[-:–]\s*([+-]?(?:\d*\.?\d+|\d+\.)(?:[eE][+-]?\d+)?)/;
    }
    /**
     * Validate the value of the slider before sending it to the back-end
     * and applying it to the other views on the page.
     */
    _validate_slide_value(x) {
        return x;
    }
}
class FloatTextModel extends FloatModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'FloatTextModel', _view_name: 'FloatTextView', disabled: false, continuous_update: false });
    }
}
class BoundedFloatTextModel extends BoundedFloatModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'BoundedFloatTextModel', _view_name: 'FloatTextView', disabled: false, continuous_update: false, step: 0.1 });
    }
}
class FloatTextView extends _widget_int__WEBPACK_IMPORTED_MODULE_1__.IntTextView {
    constructor() {
        super(...arguments);
        this._parse_value = parseFloat;
        this._default_step = 'any';
    }
    /**
     * Handle key press
     */
    handleKeypress(e) {
        // Overwrite IntTextView's handleKeypress
        // which prevents decimal points.
        e.stopPropagation();
    }
    /**
     * Handle key up
     */
    handleKeyUp(e) {
        // Overwrite IntTextView's handleKeyUp
        // which prevents decimal points.
    }
}
class FloatProgressModel extends BoundedFloatModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'FloatProgressModel', _view_name: 'ProgressView', orientation: 'horizontal', bar_style: '', style: null });
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_image.js":
/*!***************************************************!*\
  !*** ../../packages/controls/lib/widget_image.js ***!
  \***************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ImageModel: () => (/* binding */ ImageModel),
/* harmony export */   ImageView: () => (/* binding */ ImageView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


class ImageModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ImageModel', _view_name: 'ImageView', format: 'png', width: '', height: '', value: new DataView(new ArrayBuffer(0)) });
    }
}
ImageModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel.serializers), { value: {
        serialize: (value) => {
            return new DataView(value.buffer.slice(0));
        },
    } });
class ImageView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    render() {
        /**
         * Called when view is rendered.
         */
        super.render();
        this.luminoWidget.addClass('jupyter-widgets');
        this.luminoWidget.addClass('widget-image');
        this.update(); // Set defaults.
    }
    update() {
        /**
         * Update the contents of this view
         *
         * Called when the model is changed.  The model may have been
         * changed by another view or by a state update from the back-end.
         */
        let url;
        const format = this.model.get('format');
        const value = this.model.get('value');
        if (format !== 'url') {
            const blob = new Blob([value], {
                type: `image/${this.model.get('format')}`,
            });
            url = URL.createObjectURL(blob);
        }
        else {
            url = new TextDecoder('utf-8').decode(value.buffer);
        }
        // Clean up the old objectURL
        const oldurl = this.el.src;
        this.el.src = url;
        if (oldurl) {
            URL.revokeObjectURL(oldurl);
        }
        const width = this.model.get('width');
        if (width !== undefined && width.length > 0) {
            this.el.setAttribute('width', width);
        }
        else {
            this.el.removeAttribute('width');
        }
        const height = this.model.get('height');
        if (height !== undefined && height.length > 0) {
            this.el.setAttribute('height', height);
        }
        else {
            this.el.removeAttribute('height');
        }
        return super.update();
    }
    remove() {
        if (this.el.src) {
            URL.revokeObjectURL(this.el.src);
        }
        super.remove();
    }
    preinitialize() {
        // Must set this before the initialize method creates the element
        this.tagName = 'img';
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_int.js":
/*!*************************************************!*\
  !*** ../../packages/controls/lib/widget_int.js ***!
  \*************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   BaseIntSliderView: () => (/* binding */ BaseIntSliderView),
/* harmony export */   BoundedIntModel: () => (/* binding */ BoundedIntModel),
/* harmony export */   BoundedIntTextModel: () => (/* binding */ BoundedIntTextModel),
/* harmony export */   IntModel: () => (/* binding */ IntModel),
/* harmony export */   IntProgressModel: () => (/* binding */ IntProgressModel),
/* harmony export */   IntRangeSliderModel: () => (/* binding */ IntRangeSliderModel),
/* harmony export */   IntRangeSliderView: () => (/* binding */ IntRangeSliderView),
/* harmony export */   IntSliderModel: () => (/* binding */ IntSliderModel),
/* harmony export */   IntSliderView: () => (/* binding */ IntSliderView),
/* harmony export */   IntTextModel: () => (/* binding */ IntTextModel),
/* harmony export */   IntTextView: () => (/* binding */ IntTextView),
/* harmony export */   PlayModel: () => (/* binding */ PlayModel),
/* harmony export */   PlayView: () => (/* binding */ PlayView),
/* harmony export */   ProgressStyleModel: () => (/* binding */ ProgressStyleModel),
/* harmony export */   ProgressView: () => (/* binding */ ProgressView),
/* harmony export */   SliderStyleModel: () => (/* binding */ SliderStyleModel)
/* harmony export */ });
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var d3_format__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! d3-format */ "../../node_modules/d3-format/src/defaultLocale.js");
/* harmony import */ var nouislider__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! nouislider */ "../../node_modules/nouislider/dist/nouislider.js");
/* harmony import */ var nouislider__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(nouislider__WEBPACK_IMPORTED_MODULE_4__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.






class IntModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'IntModel', value: 0 });
    }
}
class BoundedIntModel extends IntModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'BoundedIntModel', max: 100, min: 0 });
    }
}
class SliderStyleModel extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'SliderStyleModel' });
    }
}
SliderStyleModel.styleProperties = Object.assign(Object.assign({}, _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel.styleProperties), { handle_color: {
        selector: '.noUi-handle',
        attribute: 'background-color',
        default: null,
    } });
class IntSliderModel extends BoundedIntModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'IntSliderModel', _view_name: 'IntSliderView', step: 1, orientation: 'horizontal', readout: true, readout_format: 'd', continuous_update: true, style: null, disabled: false });
    }
    initialize(attributes, options) {
        super.initialize(attributes, options);
        this.on('change:readout_format', this.update_readout_format, this);
        this.update_readout_format();
    }
    update_readout_format() {
        this.readout_formatter = (0,d3_format__WEBPACK_IMPORTED_MODULE_5__.format)(this.get('readout_format'));
    }
}
class IntRangeSliderModel extends IntSliderModel {
}
class BaseIntSliderView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    constructor() {
        super(...arguments);
        this._parse_value = parseInt;
    }
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-slider');
        this.el.classList.add('widget-hslider');
        // Creating noUiSlider instance and scaffolding
        this.$slider = document.createElement('div');
        this.$slider.classList.add('slider');
        // Put the slider in a container
        this.slider_container = document.createElement('div');
        this.slider_container.classList.add('slider-container');
        this.slider_container.appendChild(this.$slider);
        this.el.appendChild(this.slider_container);
        this.readout = document.createElement('div');
        this.el.appendChild(this.readout);
        this.readout.classList.add('widget-readout');
        this.readout.contentEditable = 'true';
        this.readout.style.display = 'none';
        // noUiSlider constructor and event handlers
        this.createSlider();
        // Event handlers
        this.model.on('change:orientation', this.regenSlider, this);
        this.model.on('change:max', this.updateSliderOptions, this);
        this.model.on('change:min', this.updateSliderOptions, this);
        this.model.on('change:step', this.updateSliderOptions, this);
        this.model.on('change:value', this.updateSliderValue, this);
        // Set defaults.
        this.update();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if (options === undefined || options.updated_view !== this) {
            if (this.model.get('disabled')) {
                this.readout.contentEditable = 'false';
                this.$slider.setAttribute('disabled', true);
            }
            else {
                this.readout.contentEditable = 'true';
                this.$slider.removeAttribute('disabled');
            }
            // Use the right CSS classes for vertical & horizontal sliders
            const orientation = this.model.get('orientation');
            if (orientation === 'vertical') {
                this.el.classList.remove('widget-hslider');
                this.el.classList.add('widget-vslider');
                this.el.classList.remove('widget-inline-hbox');
                this.el.classList.add('widget-inline-vbox');
            }
            else {
                this.el.classList.remove('widget-vslider');
                this.el.classList.add('widget-hslider');
                this.el.classList.remove('widget-inline-vbox');
                this.el.classList.add('widget-inline-hbox');
            }
            const readout = this.model.get('readout');
            if (readout) {
                this.readout.style.display = '';
                this.displayed.then(() => {
                    if (this.readout_overflow()) {
                        this.readout.classList.add('overflow');
                    }
                    else {
                        this.readout.classList.remove('overflow');
                    }
                });
            }
            else {
                this.readout.style.display = 'none';
            }
        }
        return super.update();
    }
    /**
     * Returns true if the readout box content overflows.
     */
    readout_overflow() {
        return this.readout.scrollWidth > this.readout.clientWidth;
    }
    events() {
        return {
            // Dictionary of events and their handlers.
            'blur [contentEditable=true]': 'handleTextChange',
            'keydown [contentEditable=true]': 'handleKeyDown',
        };
    }
    handleKeyDown(e) {
        if (e.keyCode === 13) {
            /* keyboard keycodes `enter` */
            e.preventDefault();
            this.handleTextChange();
        }
    }
    /**
     * Create a new noUiSlider object
     */
    createSlider() {
        const orientation = this.model.get('orientation');
        const behavior = this.model.get('behavior');
        nouislider__WEBPACK_IMPORTED_MODULE_4___default().create(this.$slider, {
            start: this.model.get('value'),
            connect: true,
            behaviour: behavior,
            range: {
                min: this.model.get('min'),
                max: this.model.get('max'),
            },
            step: this.model.get('step'),
            animate: false,
            orientation: orientation,
            direction: orientation === 'horizontal' ? 'ltr' : 'rtl',
            format: {
                from: (value) => Number(value),
                to: (value) => this._validate_slide_value(value),
            },
        });
        // Using noUiSlider's 'update' and 'change' events.
        // See reference: https://refreshless.com/nouislider/events-callbacks/
        this.$slider.noUiSlider.on('update', (values, handle) => {
            this.handleSliderUpdateEvent(values, handle);
        });
        this.$slider.noUiSlider.on('change', (values, handle) => {
            this.handleSliderChangeEvent(values, handle);
        });
    }
    /**
     * Recreate/Regenerate a slider object
     * noUiSlider does not support in-place mutation of the orientation
     * state. We therefore need to destroy the current instance
     * and create a new one with the new properties. This is
     * handled in a separate function and has a dedicated event
     * handler.
     */
    regenSlider(e) {
        this.$slider.noUiSlider.destroy();
        this.createSlider();
    }
    /**
     * Validate the value of the slider before sending it to the back-end
     * and applying it to the other views on the page.
     */
    _validate_slide_value(x) {
        return Math.round(x);
    }
}
class IntRangeSliderView extends BaseIntSliderView {
    constructor() {
        super(...arguments);
        // range numbers can be separated by a hyphen, colon, or an en-dash
        this._range_regex = /^\s*([+-]?\d+)\s*[-:–]\s*([+-]?\d+)/;
    }
    update(options) {
        super.update(options);
        const value = this.model.get('value');
        this.readout.textContent = this.valueToString(value);
        if (this.model.get('value') !== value) {
            this.model.set('value', value, { updated_view: this });
            this.touch();
        }
    }
    /**
     * Write value to a string
     */
    valueToString(value) {
        const format = this.model.readout_formatter;
        return value
            .map(function (v) {
            return format(v);
        })
            .join(' – ');
    }
    /**
     * Parse value from a string
     */
    stringToValue(text) {
        if (text === null) {
            return null;
        }
        // ranges can be expressed either 'val-val' or 'val:val' (+spaces)
        const match = this._range_regex.exec(text);
        if (match) {
            return [this._parse_value(match[1]), this._parse_value(match[2])];
        }
        else {
            return null;
        }
    }
    handleTextChange() {
        let value = this.stringToValue(this.readout.textContent);
        const vmin = this.model.get('min');
        const vmax = this.model.get('max');
        // reject input where NaN or lower > upper
        if (value === null ||
            isNaN(value[0]) ||
            isNaN(value[1]) ||
            value[0] > value[1]) {
            this.readout.textContent = this.valueToString(this.model.get('value'));
        }
        else {
            // clamp to range
            value = [
                Math.max(Math.min(value[0], vmax), vmin),
                Math.max(Math.min(value[1], vmax), vmin),
            ];
            if (value[0] !== this.model.get('value')[0] ||
                value[1] !== this.model.get('value')[1]) {
                this.readout.textContent = this.valueToString(value);
                this.model.set('value', value);
                this.touch();
            }
            else {
                this.readout.textContent = this.valueToString(this.model.get('value'));
            }
        }
    }
    /**
     * Called when the slider handle is released after dragging,
     * or by tapping or moving by the arrow keys.
     */
    handleSliderChangeEvent(values, handle) {
        const actual_value = values.map(this._validate_slide_value);
        this.readout.textContent = this.valueToString(actual_value);
        this.handleSliderChanged(values, handle);
    }
    /**
     * Called whilst the slider is dragged, tapped or moved by the arrow keys.
     */
    handleSliderUpdateEvent(values, handle) {
        const actual_value = values.map(this._validate_slide_value);
        this.readout.textContent = this.valueToString(actual_value);
        // Only persist the value while sliding if the continuous_update
        // trait is set to true.
        if (this.model.get('continuous_update')) {
            this.handleSliderChanged(values, handle);
        }
    }
    handleSliderChanged(values, handle) {
        const actual_value = values.map(this._validate_slide_value);
        this.model.set('value', actual_value, { updated_view: this });
        this.touch();
    }
    updateSliderOptions(e) {
        this.$slider.noUiSlider.updateOptions({
            start: this.model.get('value'),
            range: {
                min: this.model.get('min'),
                max: this.model.get('max'),
            },
            step: this.model.get('step'),
        });
    }
    updateSliderValue(model, _, options) {
        if (options.updated_view === this) {
            return;
        }
        const prev_value = this.$slider.noUiSlider.get();
        const value = this.model.get('value');
        if (prev_value[0] !== value[0] || prev_value[1] !== value[1]) {
            this.$slider.noUiSlider.set(value);
        }
    }
}
class IntSliderView extends BaseIntSliderView {
    update(options) {
        super.update(options);
        const min = this.model.get('min');
        const max = this.model.get('max');
        let value = this.model.get('value');
        if (value > max) {
            value = max;
        }
        else if (value < min) {
            value = min;
        }
        this.readout.textContent = this.valueToString(value);
        if (this.model.get('value') !== value) {
            this.model.set('value', value, { updated_view: this });
            this.touch();
        }
    }
    valueToString(value) {
        const format = this.model.readout_formatter;
        return format(value);
    }
    stringToValue(text) {
        return this._parse_value(text);
    }
    handleTextChange() {
        var _a;
        let value = this.stringToValue((_a = this.readout.textContent) !== null && _a !== void 0 ? _a : '');
        const vmin = this.model.get('min');
        const vmax = this.model.get('max');
        if (isNaN(value)) {
            this.readout.textContent = this.valueToString(this.model.get('value'));
        }
        else {
            value = Math.max(Math.min(value, vmax), vmin);
            if (value !== this.model.get('value')) {
                this.readout.textContent = this.valueToString(value);
                this.model.set('value', value);
                this.touch();
            }
            else {
                this.readout.textContent = this.valueToString(this.model.get('value'));
            }
        }
    }
    handleSliderChangeEvent(values, handle) {
        const actual_value = values.map(this._validate_slide_value);
        this.readout.textContent = this.valueToString(actual_value);
        this.handleSliderChanged(values, handle);
    }
    handleSliderUpdateEvent(values, handle) {
        const actual_value = values.map(this._validate_slide_value);
        this.readout.textContent = this.valueToString(actual_value);
        // Only persist the value while sliding if the continuous_update
        // trait is set to true.
        if (this.model.get('continuous_update')) {
            this.handleSliderChanged(values, handle);
        }
    }
    handleSliderChanged(values, handle) {
        const actual_value = this._validate_slide_value(values[handle]);
        const model_value = this.model.get('value');
        if (parseFloat(model_value) !== actual_value) {
            this.model.set('value', actual_value, { updated_view: this });
            this.touch();
        }
    }
    updateSliderOptions(e) {
        this.$slider.noUiSlider.updateOptions({
            start: this.model.get('value'),
            range: {
                min: this.model.get('min'),
                max: this.model.get('max'),
            },
            step: this.model.get('step'),
        });
    }
    updateSliderValue(model, _, options) {
        if (options.updated_view === this) {
            return;
        }
        const prev_value = this.$slider.noUiSlider.get();
        const value = this.model.get('value');
        if (prev_value !== value) {
            this.$slider.noUiSlider.set(value);
        }
    }
}
class IntTextModel extends IntModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'IntTextModel', _view_name: 'IntTextView', disabled: false, continuous_update: false });
    }
}
class BoundedIntTextModel extends BoundedIntModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'BoundedIntTextModel', _view_name: 'IntTextView', disabled: false, continuous_update: false, step: 1 });
    }
}
class IntTextView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    constructor() {
        super(...arguments);
        this._parse_value = parseInt;
        this._default_step = '1';
    }
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-text');
        this.textbox = document.createElement('input');
        this.textbox.type = 'number';
        this.textbox.required = true;
        this.textbox.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_3__.uuid)();
        this.el.appendChild(this.textbox);
        this.update(); // Set defaults.
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if (options === undefined || options.updated_view !== this) {
            const value = this.model.get('value');
            if (this._parse_value(this.textbox.value) !== value) {
                this.textbox.value = value.toString();
            }
            if (this.model.get('min') !== undefined) {
                this.textbox.min = this.model.get('min');
            }
            if (this.model.get('max') !== undefined) {
                this.textbox.max = this.model.get('max');
            }
            if (this.model.get('step') !== undefined &&
                this.model.get('step') !== null) {
                this.textbox.step = this.model.get('step');
            }
            else {
                this.textbox.step = this._default_step;
            }
            this.textbox.disabled = this.model.get('disabled');
        }
        return super.update();
    }
    events() {
        return {
            'keydown input': 'handleKeyDown',
            'keypress input': 'handleKeypress',
            'keyup input': 'handleKeyUp',
            'input input': 'handleChanging',
            'change input': 'handleChanged',
        };
    }
    /**
     * Handle key down
     *
     * Stop propagation so the event isn't sent to the application.
     */
    handleKeyDown(e) {
        e.stopPropagation();
    }
    /**
     * Handles key press
     */
    handleKeypress(e) {
        if (/[e,. ]/.test(String.fromCharCode(e.keyCode))) {
            e.preventDefault();
        }
    }
    /**
     * Handle key up
     */
    handleKeyUp(e) {
        if (e.altKey || e.ctrlKey) {
            return;
        }
        const target = e.target;
        /* remove invalid characters */
        let value = target.value;
        value = value.replace(/[e,.\s]/g, '');
        if (value.length >= 1) {
            const subvalue = value.substr(1);
            value = value[0] + subvalue.replace(/[+-]/g, '');
        }
        if (target.value !== value) {
            e.preventDefault();
            target.value = value;
        }
    }
    /**
     * Call the submit handler if continuous update is true and we are not
     * obviously incomplete.
     */
    handleChanging(e) {
        const target = e.target;
        const trimmed = target.value.trim();
        if (trimmed === '' || ['-', '-.', '.', '+.', '+'].indexOf(trimmed) >= 0) {
            // incomplete number
            return;
        }
        if (this.model.get('continuous_update')) {
            this.handleChanged(e);
        }
    }
    /**
     * Applies validated input.
     */
    handleChanged(e) {
        const target = e.target;
        let numericalValue = this._parse_value(target.value);
        // If parse failed, reset value to value stored in model.
        if (isNaN(numericalValue)) {
            target.value = this.model.get('value');
        }
        else {
            // Handle both the unbounded and bounded case by
            // checking to see if the max/min properties are defined
            let boundedValue = numericalValue;
            if (this.model.get('max') !== undefined) {
                boundedValue = Math.min(this.model.get('max'), boundedValue);
            }
            if (this.model.get('min') !== undefined) {
                boundedValue = Math.max(this.model.get('min'), boundedValue);
            }
            if (boundedValue !== numericalValue) {
                target.value = boundedValue;
                numericalValue = boundedValue;
            }
            // Apply the value if it has changed.
            if (numericalValue !== this.model.get('value')) {
                this.model.set('value', numericalValue, { updated_view: this });
                this.touch();
            }
        }
    }
}
class ProgressStyleModel extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ProgressStyleModel' });
    }
}
ProgressStyleModel.styleProperties = Object.assign(Object.assign({}, _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel.styleProperties), { bar_color: {
        selector: '.progress-bar',
        attribute: 'background-color',
        default: null,
    } });
class IntProgressModel extends BoundedIntModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'IntProgressModel', _view_name: 'ProgressView', orientation: 'horizontal', bar_style: '', style: null });
    }
}
class ProgressView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    initialize(parameters) {
        super.initialize(parameters);
        this.listenTo(this.model, 'change:bar_style', this.update_bar_style);
        this.luminoWidget.addClass('jupyter-widgets');
    }
    render() {
        super.render();
        const orientation = this.model.get('orientation');
        const className = orientation === 'horizontal' ? 'widget-hprogress' : 'widget-vprogress';
        this.el.classList.add(className);
        this.progress = document.createElement('div');
        this.progress.classList.add('progress');
        this.progress.style.position = 'relative';
        this.el.appendChild(this.progress);
        this.bar = document.createElement('div');
        this.bar.classList.add('progress-bar');
        this.bar.style.position = 'absolute';
        this.bar.style.bottom = '0px';
        this.bar.style.left = '0px';
        this.progress.appendChild(this.bar);
        // Set defaults.
        this.update();
        this.set_bar_style();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        const value = this.model.get('value');
        const max = this.model.get('max');
        const min = this.model.get('min');
        const orientation = this.model.get('orientation');
        const percent = (100.0 * (value - min)) / (max - min);
        if (orientation === 'horizontal') {
            this.el.classList.remove('widget-inline-vbox');
            this.el.classList.remove('widget-vprogress');
            this.el.classList.add('widget-inline-hbox');
            this.el.classList.add('widget-hprogress');
            this.bar.style.width = percent + '%';
            this.bar.style.height = '100%';
        }
        else {
            this.el.classList.remove('widget-inline-hbox');
            this.el.classList.remove('widget-hprogress');
            this.el.classList.add('widget-inline-vbox');
            this.el.classList.add('widget-vprogress');
            this.bar.style.width = '100%';
            this.bar.style.height = percent + '%';
        }
        return super.update();
    }
    update_bar_style() {
        this.update_mapped_classes(ProgressView.class_map, 'bar_style', this.bar);
    }
    set_bar_style() {
        this.set_mapped_classes(ProgressView.class_map, 'bar_style', this.bar);
    }
}
ProgressView.class_map = {
    success: ['progress-bar-success'],
    info: ['progress-bar-info'],
    warning: ['progress-bar-warning'],
    danger: ['progress-bar-danger'],
};
class PlayModel extends BoundedIntModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'PlayModel', _view_name: 'PlayView', repeat: false, playing: false, show_repeat: true, interval: 100, step: 1, disabled: false });
    }
    initialize(attributes, options) {
        super.initialize(attributes, options);
    }
    loop() {
        if (!this.get('playing')) {
            return;
        }
        const next_value = this.get('value') + this.get('step');
        if (next_value <= this.get('max')) {
            this.set('value', next_value);
            this.schedule_next();
        }
        else {
            if (this.get('repeat')) {
                this.set('value', this.get('min'));
                this.schedule_next();
            }
            else {
                this.pause();
            }
        }
        this.save_changes();
    }
    schedule_next() {
        this._timerId = window.setTimeout(this.loop.bind(this), this.get('interval'));
    }
    stop() {
        this.pause();
        this.set('value', this.get('min'));
        this.save_changes();
    }
    pause() {
        window.clearTimeout(this._timerId);
        this._timerId = undefined;
        this.set('playing', false);
        this.save_changes();
    }
    animate() {
        if (this._timerId !== undefined) {
            return;
        }
        if (this.get('value') === this.get('max')) {
            // if the value is at the end, reset it first, and then schedule the next
            this.set('value', this.get('min'));
            this.schedule_next();
            this.save_changes();
        }
        else {
            // otherwise directly start with the next value
            // loop will call save_changes in this case
            this.loop();
        }
        this.save_changes();
    }
    play() {
        this.set('playing', !this.get('playing'));
        this.save_changes();
    }
    repeat() {
        this.set('repeat', !this.get('repeat'));
        this.save_changes();
    }
}
class PlayView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_2__.DOMWidgetView {
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-play');
        this.playPauseButton = document.createElement('button');
        this.stopButton = document.createElement('button');
        this.repeatButton = document.createElement('button');
        this.playPauseButton.className = 'jupyter-button';
        this.stopButton.className = 'jupyter-button';
        this.repeatButton.className = 'jupyter-button';
        this.el.appendChild(this.playPauseButton); // Toggle button with playing
        this.el.appendChild(this.stopButton); // Disable if not playing
        this.el.appendChild(this.repeatButton); // Always enabled, but may be hidden
        const playIcon = document.createElement('i');
        playIcon.className = 'fa fa-play';
        this.playPauseButton.appendChild(playIcon);
        const stopIcon = document.createElement('i');
        stopIcon.className = 'fa fa-stop';
        this.stopButton.appendChild(stopIcon);
        const repeatIcon = document.createElement('i');
        repeatIcon.className = 'fa fa-retweet';
        this.repeatButton.appendChild(repeatIcon);
        this.playPauseButton.onclick = this.model.play.bind(this.model);
        this.stopButton.onclick = this.model.stop.bind(this.model);
        this.repeatButton.onclick = this.model.repeat.bind(this.model);
        this.listenTo(this.model, 'change:playing', this.onPlayingChanged);
        this.listenTo(this.model, 'change:repeat', this.updateRepeat);
        this.listenTo(this.model, 'change:show_repeat', this.updateRepeat);
        this.updatePlaying();
        this.updateRepeat();
        this.update();
    }
    update() {
        const disabled = this.model.get('disabled');
        this.playPauseButton.disabled = disabled;
        this.stopButton.disabled = disabled;
        this.repeatButton.disabled = disabled;
        this.updatePlaying();
    }
    onPlayingChanged() {
        this.updatePlaying();
        const previous = this.model.previous('playing');
        const current = this.model.get('playing');
        if (!previous && current) {
            this.model.animate();
        }
        else {
            this.model.pause();
        }
    }
    updatePlaying() {
        const playing = this.model.get('playing');
        const icon = this.playPauseButton.getElementsByTagName('i')[0];
        if (playing) {
            icon.className = 'fa fa-pause';
        }
        else {
            icon.className = 'fa fa-play';
        }
    }
    updateRepeat() {
        const repeat = this.model.get('repeat');
        this.repeatButton.style.display = this.model.get('show_repeat')
            ? this.playPauseButton.style.display
            : 'none';
        if (repeat) {
            this.repeatButton.classList.add('mod-active');
        }
        else {
            this.repeatButton.classList.remove('mod-active');
        }
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_link.js":
/*!**************************************************!*\
  !*** ../../packages/controls/lib/widget_link.js ***!
  \**************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DirectionalLinkModel: () => (/* binding */ DirectionalLinkModel),
/* harmony export */   LinkModel: () => (/* binding */ LinkModel)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


class DirectionalLinkModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { target: undefined, source: undefined, _model_name: 'DirectionalLinkModel' });
    }
    initialize(attributes, options) {
        super.initialize(attributes, options);
        this.on('change', this.updateBindings, this);
        this.updateBindings();
    }
    updateValue(sourceModel, sourceAttr, targetModel, targetAttr) {
        if (this._updating) {
            return;
        }
        this._updating = true;
        try {
            if (targetModel) {
                targetModel.set(targetAttr, sourceModel.get(sourceAttr));
                targetModel.save_changes();
            }
        }
        finally {
            this._updating = false;
        }
    }
    updateBindings() {
        this.cleanup();
        [this.sourceModel, this.sourceAttr] = this.get('source') || [null, null];
        [this.targetModel, this.targetAttr] = this.get('target') || [null, null];
        if (this.sourceModel) {
            this.listenTo(this.sourceModel, 'change:' + this.sourceAttr, () => {
                this.updateValue(this.sourceModel, this.sourceAttr, this.targetModel, this.targetAttr);
            });
            this.updateValue(this.sourceModel, this.sourceAttr, this.targetModel, this.targetAttr);
            this.listenToOnce(this.sourceModel, 'destroy', this.cleanup);
        }
        if (this.targetModel) {
            this.listenToOnce(this.targetModel, 'destroy', this.cleanup);
        }
    }
    cleanup() {
        // Stop listening to 'change' and 'destroy' events of the source and target
        if (this.sourceModel) {
            this.stopListening(this.sourceModel, 'change:' + this.sourceAttr, undefined);
            this.stopListening(this.sourceModel, 'destroy', undefined);
        }
        if (this.targetModel) {
            this.stopListening(this.targetModel, 'destroy', undefined);
        }
    }
}
DirectionalLinkModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreWidgetModel.serializers), { target: { deserialize: _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.unpack_models }, source: { deserialize: _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.unpack_models } });
class LinkModel extends DirectionalLinkModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'LinkModel' });
    }
    updateBindings() {
        super.updateBindings();
        if (this.targetModel) {
            this.listenTo(this.targetModel, 'change:' + this.targetAttr, () => {
                this.updateValue(this.targetModel, this.targetAttr, this.sourceModel, this.sourceAttr);
            });
        }
    }
    cleanup() {
        super.cleanup();
        if (this.targetModel) {
            this.stopListening(this.targetModel, 'change:' + this.targetAttr, undefined);
        }
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_selection.js":
/*!*******************************************************!*\
  !*** ../../packages/controls/lib/widget_selection.js ***!
  \*******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DropdownModel: () => (/* binding */ DropdownModel),
/* harmony export */   DropdownView: () => (/* binding */ DropdownView),
/* harmony export */   MultipleSelectionModel: () => (/* binding */ MultipleSelectionModel),
/* harmony export */   RadioButtonsModel: () => (/* binding */ RadioButtonsModel),
/* harmony export */   RadioButtonsView: () => (/* binding */ RadioButtonsView),
/* harmony export */   SelectModel: () => (/* binding */ SelectModel),
/* harmony export */   SelectMultipleModel: () => (/* binding */ SelectMultipleModel),
/* harmony export */   SelectMultipleView: () => (/* binding */ SelectMultipleView),
/* harmony export */   SelectView: () => (/* binding */ SelectView),
/* harmony export */   SelectionModel: () => (/* binding */ SelectionModel),
/* harmony export */   SelectionRangeSliderModel: () => (/* binding */ SelectionRangeSliderModel),
/* harmony export */   SelectionRangeSliderView: () => (/* binding */ SelectionRangeSliderView),
/* harmony export */   SelectionSliderModel: () => (/* binding */ SelectionSliderModel),
/* harmony export */   SelectionSliderView: () => (/* binding */ SelectionSliderView),
/* harmony export */   SelectionView: () => (/* binding */ SelectionView),
/* harmony export */   ToggleButtonsModel: () => (/* binding */ ToggleButtonsModel),
/* harmony export */   ToggleButtonsStyleModel: () => (/* binding */ ToggleButtonsStyleModel),
/* harmony export */   ToggleButtonsView: () => (/* binding */ ToggleButtonsView)
/* harmony export */ });
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var nouislider__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! nouislider */ "../../node_modules/nouislider/dist/nouislider.js");
/* harmony import */ var nouislider__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(nouislider__WEBPACK_IMPORTED_MODULE_3__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.





class SelectionModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'SelectionModel', index: '', _options_labels: [], disabled: false });
    }
}
class SelectionView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render(); // Incl. setting some defaults.
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        super.update();
        // Disable listbox if needed
        if (this.listbox) {
            this.listbox.disabled = this.model.get('disabled');
        }
        // Set tabindex
        this.updateTabindex();
        this.updateTooltip();
    }
    updateTabindex() {
        if (!this.listbox) {
            return; // we might be constructing the parent
        }
        const tabbable = this.model.get('tabbable');
        if (tabbable === true) {
            this.listbox.setAttribute('tabIndex', '0');
        }
        else if (tabbable === false) {
            this.listbox.setAttribute('tabIndex', '-1');
        }
        else if (tabbable === null) {
            this.listbox.removeAttribute('tabIndex');
        }
    }
    updateTooltip() {
        if (!this.listbox)
            return; // we might be constructing the parent
        const title = this.model.get('tooltip');
        if (!title) {
            this.listbox.removeAttribute('title');
        }
        else if (this.model.get('description').length === 0) {
            this.listbox.setAttribute('title', title);
        }
    }
}
class DropdownModel extends SelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'DropdownModel', _view_name: 'DropdownView', button_style: '' });
    }
}
// TODO: Make a Lumino dropdown control, wrapped in DropdownView. Also, fix
// bugs in keyboard handling. See
// https://github.com/jupyter-widgets/ipywidgets/issues/1055 and
// https://github.com/jupyter-widgets/ipywidgets/issues/1049
// For now, we subclass SelectView to provide DropdownView
// For the old code, see commit f68bfbc566f3a78a8f3350b438db8ed523ce3642
class DropdownView extends SelectionView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-dropdown');
        this.listbox = document.createElement('select');
        this.listbox.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.uuid)();
        this.el.appendChild(this.listbox);
        this._updateOptions();
        this.update();
    }
    /**
     * Update the contents of this view
     */
    update(options) {
        // Debounce set calls from ourselves:
        if ((options === null || options === void 0 ? void 0 : options.updated_view) !== this) {
            const optsChanged = this.model.hasChanged('_options_labels');
            if (optsChanged) {
                // Need to update options:
                this._updateOptions();
            }
        }
        // Select the correct element
        const index = this.model.get('index');
        this.listbox.selectedIndex = index === null ? -1 : index;
        return super.update();
    }
    _updateOptions() {
        this.listbox.textContent = '';
        const items = this.model.get('_options_labels');
        for (let i = 0; i < items.length; i++) {
            const item = items[i];
            const option = document.createElement('option');
            option.textContent = item.replace(/ /g, '\xa0'); // space -> &nbsp;
            option.setAttribute('data-value', encodeURIComponent(item));
            option.value = item;
            this.listbox.appendChild(option);
        }
    }
    events() {
        return {
            'change select': '_handle_change',
        };
    }
    /**
     * Handle when a new value is selected.
     */
    _handle_change() {
        this.model.set('index', this.listbox.selectedIndex === -1 ? null : this.listbox.selectedIndex, { updated_view: this });
        this.touch();
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do === 'focus') {
            this.listbox.focus();
        }
        else if (content.do === 'blur') {
            this.listbox.blur();
        }
    }
}
class SelectModel extends SelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'SelectModel', _view_name: 'SelectView', rows: 5 });
    }
}
class SelectView extends SelectionView {
    /**
     * Public constructor.
     */
    initialize(parameters) {
        super.initialize(parameters);
        // Create listbox here so that subclasses can modify it before it is populated in render()
        this.listbox = document.createElement('select');
    }
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-select');
        this.listbox.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.uuid)();
        this.el.appendChild(this.listbox);
        this._updateOptions();
        this.update();
        this.updateSelection();
    }
    /**
     * Update the contents of this view
     */
    update(options) {
        // Don't update options/index on set calls from ourselves:
        if ((options === null || options === void 0 ? void 0 : options.updated_view) !== this) {
            const optsChange = this.model.hasChanged('_options_labels');
            const idxChange = this.model.hasChanged('index');
            if (optsChange || idxChange) {
                // Stash the index to guard against change events
                const idx = this.model.get('index');
                if (optsChange) {
                    this._updateOptions();
                }
                this.updateSelection(idx);
            }
        }
        super.update();
        let rows = this.model.get('rows');
        if (rows === null) {
            rows = '';
        }
        this.listbox.setAttribute('size', rows);
    }
    updateSelection(index) {
        index = index || this.model.get('index');
        this.listbox.selectedIndex = index === null ? -1 : index;
    }
    _updateOptions() {
        this.listbox.textContent = '';
        const items = this.model.get('_options_labels');
        for (let i = 0; i < items.length; i++) {
            const item = items[i];
            const option = document.createElement('option');
            option.textContent = item.replace(/ /g, '\xa0'); // space -> &nbsp;
            option.setAttribute('data-value', encodeURIComponent(item));
            option.value = item;
            this.listbox.appendChild(option);
        }
    }
    events() {
        return {
            'change select': '_handle_change',
        };
    }
    /**
     * Handle when a new value is selected.
     */
    _handle_change() {
        this.model.set('index', this.listbox.selectedIndex, { updated_view: this });
        this.touch();
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do == 'focus') {
            this.listbox.focus();
        }
        else if (content.do == 'blur') {
            this.listbox.blur();
        }
    }
}
class RadioButtonsModel extends SelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'RadioButtonsModel', _view_name: 'RadioButtonsView', tooltips: [], icons: [], button_style: '', orientation: 'vertical' });
    }
}
class RadioButtonsView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-radio');
        this.container = document.createElement('div');
        this.el.appendChild(this.container);
        this.container.classList.add('widget-radio-box');
        this.update();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if (this.model.get('orientation') === 'vertical') {
            this.container.classList.remove('widget-radio-box-horizontal');
            this.container.classList.add('widget-radio-box-vertical');
        }
        else {
            this.container.classList.remove('widget-radio-box-vertical');
            this.container.classList.add('widget-radio-box-horizontal');
        }
        const items = this.model.get('_options_labels');
        const radios = Array.from(this.container.querySelectorAll('input[type="radio"]')).map((x) => x.value);
        let stale = items.length !== radios.length;
        if (!stale) {
            for (let i = 0, len = items.length; i < len; ++i) {
                if (radios[i] !== items[i]) {
                    stale = true;
                    break;
                }
            }
        }
        if (stale && (options === undefined || options.updated_view !== this)) {
            // Add items to the DOM.
            this.container.textContent = '';
            items.forEach((item, index) => {
                const label = document.createElement('label');
                label.textContent = item;
                this.container.appendChild(label);
                const radio = document.createElement('input');
                radio.setAttribute('type', 'radio');
                radio.value = index.toString();
                radio.setAttribute('data-value', encodeURIComponent(item));
                label.appendChild(radio);
            });
        }
        items.forEach((item, index) => {
            const item_query = 'input[data-value="' + encodeURIComponent(item) + '"]';
            const radio = this.container.querySelectorAll(item_query);
            if (radio.length > 0) {
                const radio_el = radio[0];
                radio_el.checked = this.model.get('index') === index;
                radio_el.disabled = this.model.get('disabled');
            }
        });
        // Schedule adjustPadding asynchronously to
        // allow dom elements to be created properly
        setTimeout(this.adjustPadding, 0, this);
        return super.update(options);
    }
    /**
     * Adjust Padding to Multiple of Line Height
     *
     * Adjust margins so that the overall height
     * is a multiple of a single line height.
     *
     * This widget needs it because radio options
     * are spaced tighter than individual widgets
     * yet we would like the full widget line up properly
     * when displayed side-by-side with other widgets.
     */
    adjustPadding(e) {
        // Vertical margins on a widget
        const elStyles = window.getComputedStyle(e.el);
        const margins = parseInt(elStyles.marginTop, 10) + parseInt(elStyles.marginBottom, 10);
        // Total spaces taken by a single-line widget
        const lineHeight = e.label.offsetHeight + margins;
        // Current adjustment value on this widget
        const cStyles = window.getComputedStyle(e.container);
        const containerMargin = parseInt(cStyles.marginBottom, 10);
        // How far we are off from a multiple of single windget lines
        const diff = (e.el.offsetHeight + margins - containerMargin) % lineHeight;
        // Apply the new adjustment
        const extraMargin = diff === 0 ? 0 : lineHeight - diff;
        e.container.style.marginBottom = extraMargin + 'px';
    }
    events() {
        return {
            'click input[type="radio"]': '_handle_click',
        };
    }
    /**
     * Handle when a value is clicked.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    _handle_click(event) {
        const target = event.target;
        this.model.set('index', parseInt(target.value, 10), { updated_view: this });
        this.touch();
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do == 'focus') {
            const firstItem = this.container.firstElementChild;
            firstItem.focus();
        }
        else if (content.do == 'blur') {
            for (let i = 0; i < this.container.children.length; i++) {
                const item = this.container.children[i];
                item.blur();
            }
        }
    }
}
class ToggleButtonsStyleModel extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ToggleButtonsStyleModel' });
    }
}
ToggleButtonsStyleModel.styleProperties = Object.assign(Object.assign({}, _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel.styleProperties), { button_width: {
        selector: '.widget-toggle-button',
        attribute: 'width',
        default: null,
    }, font_weight: {
        selector: '.widget-toggle-button',
        attribute: 'font-weight',
        default: '',
    } });
class ToggleButtonsModel extends SelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ToggleButtonsModel', _view_name: 'ToggleButtonsView' });
    }
}
class ToggleButtonsView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    initialize(options) {
        this._css_state = {};
        super.initialize(options);
        this.listenTo(this.model, 'change:button_style', this.update_button_style);
    }
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-toggle-buttons');
        this.buttongroup = document.createElement('div');
        this.el.appendChild(this.buttongroup);
        this.update();
        this.set_button_style();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        const items = this.model.get('_options_labels');
        const icons = this.model.get('icons') || [];
        const previous_icons = this.model.previous('icons') || [];
        const previous_bstyle = ToggleButtonsView.classMap[this.model.previous('button_style')] || '';
        const tooltips = this.model.get('tooltips') || [];
        const disabled = this.model.get('disabled');
        const buttons = this.buttongroup.querySelectorAll('button');
        const values = Array.from(buttons).map((x) => x.value);
        let stale = false;
        for (let i = 0, len = items.length; i < len; ++i) {
            if (values[i] !== items[i] || icons[i] !== previous_icons[i]) {
                stale = true;
                break;
            }
        }
        if (stale && (options === undefined || options.updated_view !== this)) {
            // Add items to the DOM.
            this.buttongroup.textContent = '';
            items.forEach((item, index) => {
                let item_html;
                const empty = item.trim().length === 0 &&
                    (!icons[index] || icons[index].trim().length === 0);
                if (empty) {
                    item_html = '&nbsp;';
                }
                else {
                    item_html = _utils__WEBPACK_IMPORTED_MODULE_2__.escape_html(item);
                }
                const icon = document.createElement('i');
                const button = document.createElement('button');
                if (icons[index]) {
                    icon.className = 'fa fa-' + icons[index];
                }
                button.setAttribute('type', 'button');
                button.className = 'widget-toggle-button jupyter-button';
                if (previous_bstyle) {
                    button.classList.add(previous_bstyle);
                }
                button.innerHTML = item_html;
                button.setAttribute('data-value', encodeURIComponent(item));
                button.setAttribute('value', index.toString());
                button.appendChild(icon);
                button.disabled = disabled;
                if (tooltips[index]) {
                    button.setAttribute('title', tooltips[index]);
                }
                this.update_style_traits(button);
                this.buttongroup.appendChild(button);
            });
        }
        // Select active button.
        items.forEach((item, index) => {
            const item_query = '[data-value="' + encodeURIComponent(item) + '"]';
            const button = this.buttongroup.querySelector(item_query);
            if (this.model.get('index') === index) {
                button.classList.add('mod-active');
            }
            else {
                button.classList.remove('mod-active');
            }
        });
        this.stylePromise.then(function (style) {
            if (style) {
                style.style();
            }
        });
        return super.update(options);
    }
    update_style_traits(button) {
        for (const name in this._css_state) {
            if (Object.prototype.hasOwnProperty.call(this._css_state, 'name')) {
                if (name === 'margin') {
                    this.buttongroup.style[name] = this._css_state[name];
                }
                else if (name !== 'width') {
                    if (button) {
                        button.style[name] = this._css_state[name];
                    }
                    else {
                        const buttons = this.buttongroup.querySelectorAll('button');
                        if (buttons.length) {
                            buttons[0].style[name] = this._css_state[name];
                        }
                    }
                }
            }
        }
    }
    update_button_style() {
        const buttons = this.buttongroup.querySelectorAll('button');
        for (let i = 0; i < buttons.length; i++) {
            this.update_mapped_classes(ToggleButtonsView.classMap, 'button_style', buttons[i]);
        }
    }
    set_button_style() {
        const buttons = this.buttongroup.querySelectorAll('button');
        for (let i = 0; i < buttons.length; i++) {
            this.set_mapped_classes(ToggleButtonsView.classMap, 'button_style', buttons[i]);
        }
    }
    events() {
        return {
            'click button': '_handle_click',
        };
    }
    /**
     * Handle when a value is clicked.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    _handle_click(event) {
        const target = event.target;
        this.model.set('index', parseInt(target.value, 10), { updated_view: this });
        this.touch();
        // We also send a clicked event, since the value is only set if it changed.
        // See https://github.com/jupyter-widgets/ipywidgets/issues/763
        this.send({ event: 'click' });
    }
}
(function (ToggleButtonsView) {
    ToggleButtonsView.classMap = {
        primary: ['mod-primary'],
        success: ['mod-success'],
        info: ['mod-info'],
        warning: ['mod-warning'],
        danger: ['mod-danger'],
    };
})(ToggleButtonsView || (ToggleButtonsView = {}));
class SelectionSliderModel extends SelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'SelectionSliderModel', _view_name: 'SelectionSliderView', orientation: 'horizontal', readout: true, continuous_update: true });
    }
}
class SelectionSliderView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-hslider');
        this.el.classList.add('widget-slider');
        // Creating noUiSlider instance and scaffolding
        this.$slider = document.createElement('div');
        this.$slider.classList.add('slider');
        // Put the slider in a container
        this.slider_container = document.createElement('div');
        this.slider_container.classList.add('slider-container');
        this.slider_container.appendChild(this.$slider);
        this.el.appendChild(this.slider_container);
        this.readout = document.createElement('div');
        this.el.appendChild(this.readout);
        this.readout.classList.add('widget-readout');
        this.readout.style.display = 'none';
        // noUiSlider constructor and event handlers
        this.createSlider();
        // Event handlers
        this.model.on('change:orientation', this.regenSlider, this);
        this.model.on('change:index', this.updateSliderValue, this);
        // Set defaults.
        this.update();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if ((options === null || options === void 0 ? void 0 : options.updated_view) !== this) {
            this.updateSliderOptions(this.model);
            const orientation = this.model.get('orientation');
            const disabled = this.model.get('disabled');
            if (disabled) {
                this.readout.contentEditable = 'false';
                this.$slider.setAttribute('disabled', true);
            }
            else {
                this.readout.contentEditable = 'true';
                this.$slider.removeAttribute('disabled');
            }
            // Use the right CSS classes for vertical & horizontal sliders
            if (orientation === 'vertical') {
                this.el.classList.remove('widget-hslider');
                this.el.classList.remove('widget-inline-hbox');
                this.el.classList.add('widget-vslider');
                this.el.classList.add('widget-inline-vbox');
            }
            else {
                this.el.classList.remove('widget-vslider');
                this.el.classList.remove('widget-inline-vbox');
                this.el.classList.add('widget-hslider');
                this.el.classList.add('widget-inline-hbox');
            }
            const readout = this.model.get('readout');
            if (readout) {
                // this.$readout.show();
                this.readout.style.display = '';
            }
            else {
                // this.$readout.hide();
                this.readout.style.display = 'none';
            }
            this.updateSelection();
        }
        return super.update(options);
    }
    regenSlider(e) {
        this.$slider.noUiSlider.destroy();
        this.createSlider();
    }
    createSlider() {
        const labels = this.model.get('_options_labels');
        const min = 0;
        const max = labels.length - 1;
        const orientation = this.model.get('orientation');
        const behavior = this.model.get('behavior');
        nouislider__WEBPACK_IMPORTED_MODULE_3___default().create(this.$slider, {
            start: this.model.get('index'),
            connect: true,
            behaviour: behavior,
            range: {
                min: min,
                max: max,
            },
            step: 1,
            animate: false,
            orientation: orientation,
            direction: orientation === 'horizontal' ? 'ltr' : 'rtl',
            format: {
                from: (value) => Number(value),
                to: (value) => Math.round(value),
            },
        });
        // Using noUiSlider's 'update' and 'change' events.
        // See reference: https://refreshless.com/nouislider/events-callbacks/
        this.$slider.noUiSlider.on('update', (values, handle) => {
            this.handleSliderUpdateEvent(values, handle);
        });
        this.$slider.noUiSlider.on('change', (values, handle) => {
            this.handleSliderChangeEvent(values, handle);
        });
    }
    events() {
        return {
            slide: 'handleSliderChange',
            slidestop: 'handleSliderChanged',
        };
    }
    updateSelection() {
        const index = this.model.get('index');
        this.updateReadout(index);
    }
    updateReadout(index) {
        const value = this.model.get('_options_labels')[index];
        this.readout.textContent = value;
    }
    /**
     * Called whilst the slider is dragged, tapped or moved by the arrow keys.
     */
    handleSliderUpdateEvent(values, handle) {
        const index = values[0];
        this.updateReadout(index);
        // Only persist the value while sliding if the continuous_update
        // trait is set to true.
        if (this.model.get('continuous_update')) {
            this.handleSliderChanged(values, handle);
        }
    }
    /**
     * Called when the slider handle is released after dragging,
     * or by tapping or moving by the arrow keys.
     */
    handleSliderChangeEvent(values, handle) {
        const index = values[0];
        this.updateReadout(index);
        this.handleSliderChanged(values, handle);
    }
    /**
     * Called when the slider value has changed.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    handleSliderChanged(values, handle) {
        const index = values[0];
        this.updateReadout(index);
        this.model.set('index', index, { updated_view: this });
        this.touch();
    }
    updateSliderOptions(e) {
        const labels = this.model.get('_options_labels');
        const min = 0;
        const max = labels.length - 1;
        this.$slider.noUiSlider.updateOptions({
            start: this.model.get('index'),
            range: {
                min: min,
                max: max,
            },
            step: 1,
        });
    }
    updateSliderValue(model, _, options) {
        if (options.updated_view === this) {
            return;
        }
        const prev_index = this.$slider.noUiSlider.get();
        const index = this.model.get('index');
        if (prev_index !== index) {
            this.$slider.noUiSlider.set(index);
        }
    }
}
class MultipleSelectionModel extends SelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'MultipleSelectionModel' });
    }
}
class SelectMultipleModel extends MultipleSelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'SelectMultipleModel', _view_name: 'SelectMultipleView', rows: null });
    }
}
class SelectMultipleView extends SelectView {
    /**
     * Public constructor.
     */
    initialize(parameters) {
        super.initialize(parameters);
        this.listbox.multiple = true;
    }
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-select-multiple');
    }
    updateSelection() {
        const selected = this.model.get('index') || [];
        const listboxOptions = this.listbox.options;
        // Clear the selection
        this.listbox.selectedIndex = -1;
        // Select the appropriate options
        selected.forEach((i) => {
            listboxOptions[i].selected = true;
        });
    }
    /**
     * Handle when a new value is selected.
     */
    _handle_change() {
        const index = Array.prototype.map.call(this.listbox.selectedOptions || [], function (option) {
            return option.index;
        });
        this.model.set('index', index, { updated_view: this });
        this.touch();
    }
}
class SelectionRangeSliderModel extends MultipleSelectionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'SelectionSliderModel', _view_name: 'SelectionSliderView', orientation: 'horizontal', readout: true, continuous_update: true });
    }
}
class SelectionRangeSliderView extends SelectionSliderView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
    }
    updateSelection(index) {
        index = index || this.model.get('index');
        this.updateReadout(index);
    }
    updateReadout(index) {
        const labels = this.model.get('_options_labels');
        const minValue = labels[index[0]];
        const maxValue = labels[index[1]];
        this.readout.textContent = `${minValue}-${maxValue}`;
    }
    /**
     * Called when the slider value is changing.
     */
    handleSliderUpdateEvent(values, handle) {
        const intValues = values.map(Math.trunc);
        this.updateReadout(intValues);
        // Only persist the value while sliding if the continuous_update
        // trait is set to true.
        if (this.model.get('continuous_update')) {
            this.handleSliderChanged(values, handle);
        }
    }
    /**
     * Called when the slider value has changed.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    handleSliderChanged(values, handle) {
        const intValues = values.map(Math.round);
        this.updateReadout(intValues);
        // set index to a snapshot of the values passed by the slider
        this.model.set('index', intValues.slice(), { updated_view: this });
        this.touch();
    }
    updateSliderValue(model, _, options) {
        if (options.updated_view === this) {
            return;
        }
        // Rounding values to avoid floating point precision error for the if statement below
        const prev_index = this.$slider.noUiSlider.get().map(Math.round);
        const index = this.model.get('index').map(Math.round);
        if (prev_index[0] !== index[0] || prev_index[1] !== index[1]) {
            this.$slider.noUiSlider.set(index);
        }
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_selectioncontainer.js":
/*!****************************************************************!*\
  !*** ../../packages/controls/lib/widget_selectioncontainer.js ***!
  \****************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   AccordionModel: () => (/* binding */ AccordionModel),
/* harmony export */   AccordionView: () => (/* binding */ AccordionView),
/* harmony export */   JupyterLuminoAccordionWidget: () => (/* binding */ JupyterLuminoAccordionWidget),
/* harmony export */   JupyterLuminoTabPanelWidget: () => (/* binding */ JupyterLuminoTabPanelWidget),
/* harmony export */   SelectionContainerModel: () => (/* binding */ SelectionContainerModel),
/* harmony export */   StackModel: () => (/* binding */ StackModel),
/* harmony export */   StackView: () => (/* binding */ StackView),
/* harmony export */   TabModel: () => (/* binding */ TabModel),
/* harmony export */   TabView: () => (/* binding */ TabView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_box__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_box */ "../../packages/controls/lib/widget_box.js");
/* harmony import */ var _lumino_tabpanel__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./lumino/tabpanel */ "../../packages/controls/lib/lumino/tabpanel.js");
/* harmony import */ var _lumino_accordion__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./lumino/accordion */ "../../packages/controls/lib/lumino/accordion.js");
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @lumino/widgets */ "webpack/sharing/consume/default/@lumino/widgets");
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(_lumino_widgets__WEBPACK_IMPORTED_MODULE_4__);
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @lumino/algorithm */ "webpack/sharing/consume/default/@lumino/algorithm");
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_5___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_5__);
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! @lumino/messaging */ "webpack/sharing/consume/default/@lumino/messaging");
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_6___default = /*#__PURE__*/__webpack_require__.n(_lumino_messaging__WEBPACK_IMPORTED_MODULE_6__);
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! jquery */ "webpack/sharing/consume/default/jquery/jquery?123b");
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_7___default = /*#__PURE__*/__webpack_require__.n(jquery__WEBPACK_IMPORTED_MODULE_7__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.








class SelectionContainerModel extends _widget_box__WEBPACK_IMPORTED_MODULE_1__.BoxModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'SelectionContainerModel', selected_index: null, titles: [] });
    }
}
class AccordionModel extends SelectionContainerModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'AccordionModel', _view_name: 'AccordionView' });
    }
}
// We implement our own tab widget since Phoshpor's TabPanel uses an absolute
// positioning BoxLayout, but we want a more an html/css-based Panel layout.
class JupyterLuminoAccordionWidget extends _lumino_accordion__WEBPACK_IMPORTED_MODULE_3__.Accordion {
    constructor(options) {
        const view = options.view;
        delete options.view;
        super(options);
        this._view = view;
    }
    /**
     * Process the Lumino message.
     *
     * Any custom Lumino widget used inside a Jupyter widget should override
     * the processMessage function like this.
     */
    processMessage(msg) {
        var _a;
        super.processMessage(msg);
        (_a = this._view) === null || _a === void 0 ? void 0 : _a.processLuminoMessage(msg);
    }
    /**
     * Dispose the widget.
     *
     * This causes the view to be destroyed as well with 'remove'
     */
    dispose() {
        if (this.isDisposed) {
            return;
        }
        super.dispose();
        this._view.remove();
        this._view = null;
    }
}
class AccordionView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    _createElement(tagName) {
        this.luminoWidget = new JupyterLuminoAccordionWidget({ view: this });
        return this.luminoWidget.node;
    }
    _setElement(el) {
        if (this.el || el !== this.luminoWidget.node) {
            // Accordions don't allow setting the element beyond the initial creation.
            throw new Error('Cannot reset the DOM element.');
        }
        this.el = this.luminoWidget.node;
        this.$el = jquery__WEBPACK_IMPORTED_MODULE_7___default()(this.luminoWidget.node);
    }
    initialize(parameters) {
        super.initialize(parameters);
        this.children_views = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.ViewList(this.add_child_view, this.remove_child_view, this);
        this.listenTo(this.model, 'change:children', () => this.updateChildren());
        this.listenTo(this.model, 'change:selected_index', () => this.update_selected_index());
        this.listenTo(this.model, 'change:titles', () => this.update_titles());
    }
    /**
     * Called when view is rendered.
     */
    render() {
        var _a;
        super.render();
        const accordion = this.luminoWidget;
        accordion.addClass('jupyter-widgets');
        accordion.addClass('widget-accordion');
        accordion.addClass('widget-container');
        accordion.selection.selectionChanged.connect((sender) => {
            if (!this.updatingChildren) {
                this.model.set('selected_index', accordion.selection.index);
                this.touch();
            }
        });
        (_a = this.children_views) === null || _a === void 0 ? void 0 : _a.update(this.model.get('children'));
        this.update_titles();
        this.update_selected_index();
    }
    /**
     * Update children
     */
    updateChildren() {
        var _a;
        // While we are updating, the index may not be valid, so deselect the
        // tabs before updating so we don't get spurious changes in the index,
        // which would then set off another sync cycle.
        this.updatingChildren = true;
        this.luminoWidget.selection.index = null;
        (_a = this.children_views) === null || _a === void 0 ? void 0 : _a.update(this.model.get('children'));
        this.update_selected_index();
        this.updatingChildren = false;
    }
    /**
     * Set header titles
     */
    update_titles() {
        const collapsed = this.luminoWidget.collapseWidgets;
        const titles = this.model.get('titles');
        for (let i = 0; i < collapsed.length; i++) {
            if (titles[i] !== void 0) {
                collapsed[i].widget.title.label = titles[i];
            }
        }
    }
    /**
     * Make the rendering and selected index consistent.
     */
    update_selected_index() {
        this.luminoWidget.selection.index = this.model.get('selected_index');
    }
    /**
     * Called when a child is removed from children list.
     */
    remove_child_view(view) {
        this.luminoWidget.removeWidget(view.luminoWidget);
        view.remove();
    }
    /**
     * Called when a child is added to children list.
     */
    add_child_view(model, index) {
        // Placeholder widget to keep our position in the tab panel while we create the view.
        const accordion = this.luminoWidget;
        const placeholder = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_4__.Widget();
        placeholder.title.label = this.model.get('titles')[index] || '';
        accordion.addWidget(placeholder);
        return this.create_child_view(model)
            .then((view) => {
            const widget = view.luminoWidget;
            widget.title.label = placeholder.title.label;
            const collapse = accordion.collapseWidgets[accordion.indexOf(placeholder)];
            collapse.widget = widget;
            placeholder.dispose();
            return view;
        })
            .catch((0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.reject)('Could not add child view to box', true));
    }
    remove() {
        this.children_views = null;
        super.remove();
    }
}
class TabModel extends SelectionContainerModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'TabModel', _view_name: 'TabView' });
    }
}
// We implement our own tab widget since Phoshpor's TabPanel uses an absolute
// positioning BoxLayout, but we want a more an html/css-based Panel layout.
class JupyterLuminoTabPanelWidget extends _lumino_tabpanel__WEBPACK_IMPORTED_MODULE_2__.TabPanel {
    constructor(options) {
        const view = options.view;
        delete options.view;
        super(options);
        this._view = view;
        // We want the view's messages to be the messages the tabContents panel
        // gets.
        _lumino_messaging__WEBPACK_IMPORTED_MODULE_6__.MessageLoop.installMessageHook(this.tabContents, (handler, msg) => {
            // There may be times when we want the view's handler to be called
            // *after* the message has been processed by the widget, in which
            // case we'll need to revisit using a message hook.
            this._view.processLuminoMessage(msg);
            return true;
        });
    }
    /**
     * Dispose the widget.
     *
     * This causes the view to be destroyed as well with 'remove'
     */
    dispose() {
        if (this.isDisposed) {
            return;
        }
        super.dispose();
        this._view.remove();
        this._view = null;
    }
}
class TabView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    constructor() {
        super(...arguments);
        this.updatingTabs = false;
    }
    _createElement(tagName) {
        this.luminoWidget = new JupyterLuminoTabPanelWidget({
            view: this,
        });
        return this.luminoWidget.node;
    }
    _setElement(el) {
        if (this.el || el !== this.luminoWidget.node) {
            // TabViews don't allow setting the element beyond the initial creation.
            throw new Error('Cannot reset the DOM element.');
        }
        this.el = this.luminoWidget.node;
        this.$el = jquery__WEBPACK_IMPORTED_MODULE_7___default()(this.luminoWidget.node);
    }
    /**
     * Public constructor.
     */
    initialize(parameters) {
        super.initialize(parameters);
        this.childrenViews = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.ViewList(this.addChildView, (view) => {
            view.remove();
        }, this);
        this.listenTo(this.model, 'change:children', () => this.updateTabs());
        this.listenTo(this.model, 'change:titles', () => this.updateTitles());
    }
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        const tabs = this.luminoWidget;
        tabs.addClass('jupyter-widgets');
        tabs.addClass('widget-container');
        tabs.addClass('jupyter-widget-tab');
        tabs.addClass('widget-tab');
        tabs.tabsMovable = true;
        tabs.tabBar.insertBehavior = 'none'; // needed for insert behavior, see below.
        tabs.tabBar.currentChanged.connect(this._onTabChanged, this);
        tabs.tabBar.tabMoved.connect(this._onTabMoved, this);
        tabs.tabBar.addClass('widget-tab-bar');
        tabs.tabContents.addClass('widget-tab-contents');
        // TODO: expose this option in python
        tabs.tabBar.tabsMovable = false;
        this.updateTabs();
        this.update();
    }
    /**
     * Render tab views based on the current model's children.
     */
    updateTabs() {
        var _a;
        // While we are updating, the index may not be valid, so deselect the
        // tabs before updating so we don't get spurious changes in the index,
        // which would then set off another sync cycle.
        this.updatingTabs = true;
        this.luminoWidget.currentIndex = null;
        (_a = this.childrenViews) === null || _a === void 0 ? void 0 : _a.update(this.model.get('children'));
        this.luminoWidget.currentIndex = this.model.get('selected_index');
        this.updatingTabs = false;
    }
    /**
     * Called when a child is added to children list.
     */
    addChildView(model, index) {
        // Placeholder widget to keep our position in the tab panel while we create the view.
        const label = this.model.get('titles')[index] || '';
        const tabs = this.luminoWidget;
        const placeholder = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_4__.Widget();
        placeholder.title.label = label;
        tabs.addWidget(placeholder);
        return this.create_child_view(model)
            .then((view) => {
            const widget = view.luminoWidget;
            widget.title.label = placeholder.title.label;
            widget.title.closable = false;
            const i = _lumino_algorithm__WEBPACK_IMPORTED_MODULE_5__.ArrayExt.firstIndexOf(tabs.widgets, placeholder);
            // insert after placeholder so that if placeholder is selected, the
            // real widget will be selected now (this depends on the tab bar
            // insert behavior)
            tabs.insertWidget(i + 1, widget);
            placeholder.dispose();
            return view;
        })
            .catch((0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.reject)('Could not add child view to box', true));
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        // Update the selected index in the overall update method because it
        // should be run after the tabs have been updated. Otherwise the
        // selected index may not be a valid tab in the tab bar.
        this.updateSelectedIndex();
        return super.update();
    }
    /**
     * Updates the tab page titles.
     */
    updateTitles() {
        const titles = this.model.get('titles') || [];
        (0,_lumino_algorithm__WEBPACK_IMPORTED_MODULE_5__.each)(this.luminoWidget.widgets, (widget, i) => {
            widget.title.label = titles[i] || '';
        });
    }
    /**
     * Updates the selected index.
     */
    updateSelectedIndex() {
        this.luminoWidget.currentIndex = this.model.get('selected_index');
    }
    remove() {
        this.childrenViews = null;
        super.remove();
    }
    _onTabChanged(sender, args) {
        if (!this.updatingTabs) {
            const i = args.currentIndex;
            this.model.set('selected_index', i === -1 ? null : i);
            this.touch();
        }
    }
    /**
     * Handle the `tabMoved` signal from the tab bar.
     */
    _onTabMoved(sender, args) {
        const children = this.model.get('children').slice();
        _lumino_algorithm__WEBPACK_IMPORTED_MODULE_5__.ArrayExt.move(children, args.fromIndex, args.toIndex);
        this.model.set('children', children);
        this.touch();
    }
}
class StackModel extends SelectionContainerModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'StackModel', _view_name: 'StackView' });
    }
}
class StackView extends _widget_box__WEBPACK_IMPORTED_MODULE_1__.BoxView {
    initialize(parameters) {
        super.initialize(parameters);
        this.listenTo(this.model, 'change:selected_index', this.update_children);
    }
    update_children() {
        var _a;
        let child;
        if (this.model.get('selected_index') === null) {
            child = [];
        }
        else {
            child = [this.model.get('children')[this.model.get('selected_index')]];
        }
        (_a = this.children_views) === null || _a === void 0 ? void 0 : _a.update(child).then((views) => {
            // Notify all children that their sizes may have changed.
            views.forEach((view) => {
                _lumino_messaging__WEBPACK_IMPORTED_MODULE_6__.MessageLoop.postMessage(view.luminoWidget, _lumino_widgets__WEBPACK_IMPORTED_MODULE_4__.Widget.ResizeMessage.UnknownSize);
            });
        });
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_string.js":
/*!****************************************************!*\
  !*** ../../packages/controls/lib/widget_string.js ***!
  \****************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ComboboxModel: () => (/* binding */ ComboboxModel),
/* harmony export */   ComboboxView: () => (/* binding */ ComboboxView),
/* harmony export */   HTMLMathModel: () => (/* binding */ HTMLMathModel),
/* harmony export */   HTMLMathStyleModel: () => (/* binding */ HTMLMathStyleModel),
/* harmony export */   HTMLMathView: () => (/* binding */ HTMLMathView),
/* harmony export */   HTMLModel: () => (/* binding */ HTMLModel),
/* harmony export */   HTMLStyleModel: () => (/* binding */ HTMLStyleModel),
/* harmony export */   HTMLView: () => (/* binding */ HTMLView),
/* harmony export */   LabelModel: () => (/* binding */ LabelModel),
/* harmony export */   LabelStyleModel: () => (/* binding */ LabelStyleModel),
/* harmony export */   LabelView: () => (/* binding */ LabelView),
/* harmony export */   PasswordModel: () => (/* binding */ PasswordModel),
/* harmony export */   PasswordView: () => (/* binding */ PasswordView),
/* harmony export */   StringModel: () => (/* binding */ StringModel),
/* harmony export */   StringView: () => (/* binding */ StringView),
/* harmony export */   TextModel: () => (/* binding */ TextModel),
/* harmony export */   TextStyleModel: () => (/* binding */ TextStyleModel),
/* harmony export */   TextView: () => (/* binding */ TextView),
/* harmony export */   TextareaModel: () => (/* binding */ TextareaModel),
/* harmony export */   TextareaView: () => (/* binding */ TextareaView)
/* harmony export */ });
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./version */ "../../packages/controls/lib/version.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.




/**
 * Class name for a combobox with an invalid value.
 */
const INVALID_VALUE_CLASS = 'jpwidgets-invalidComboValue';
class StringStyleModel extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'StringStyleModel', _model_module: '@jupyter-widgets/controls', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_3__.JUPYTER_CONTROLS_VERSION });
    }
}
StringStyleModel.styleProperties = Object.assign(Object.assign({}, _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel.styleProperties), { background: {
        selector: '',
        attribute: 'background',
        default: null,
    }, font_size: {
        selector: '',
        attribute: 'font-size',
        default: '',
    }, text_color: {
        selector: '',
        attribute: 'color',
        default: '',
    } });
class HTMLStyleModel extends StringStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'HTMLStyleModel', _model_module: '@jupyter-widgets/controls', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_3__.JUPYTER_CONTROLS_VERSION });
    }
}
HTMLStyleModel.styleProperties = Object.assign({}, StringStyleModel.styleProperties);
class HTMLMathStyleModel extends StringStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'HTMLMathStyleModel', _model_module: '@jupyter-widgets/controls', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_3__.JUPYTER_CONTROLS_VERSION });
    }
}
HTMLMathStyleModel.styleProperties = Object.assign({}, StringStyleModel.styleProperties);
class LabelStyleModel extends StringStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'LabelStyleModel', _model_module: '@jupyter-widgets/controls', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_3__.JUPYTER_CONTROLS_VERSION });
    }
}
LabelStyleModel.styleProperties = Object.assign(Object.assign({}, StringStyleModel.styleProperties), { font_family: {
        selector: '',
        attribute: 'font-family',
        default: '',
    }, font_style: {
        selector: '',
        attribute: 'font-style',
        default: '',
    }, font_variant: {
        selector: '',
        attribute: 'font-variant',
        default: '',
    }, font_weight: {
        selector: '',
        attribute: 'font-weight',
        default: '',
    }, text_decoration: {
        selector: '',
        attribute: 'text-decoration',
        default: '',
    } });
class TextStyleModel extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'TextStyleModel', _model_module: '@jupyter-widgets/controls', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_3__.JUPYTER_CONTROLS_VERSION });
    }
}
TextStyleModel.styleProperties = Object.assign(Object.assign({}, _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionStyleModel.styleProperties), { background: {
        selector: '.widget-input',
        attribute: 'background',
        default: null,
    }, font_size: {
        selector: '.widget-input',
        attribute: 'font-size',
        default: '',
    }, text_color: {
        selector: '.widget-input',
        attribute: 'color',
        default: '',
    } });
class StringModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { value: '', disabled: false, placeholder: '\u200b', _model_name: 'StringModel' });
    }
}
class StringView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render(); // Incl. setting some defaults.
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
    }
}
class HTMLModel extends StringModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'HTMLView', _model_name: 'HTMLModel' });
    }
}
class HTMLView extends StringView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-html');
        this.content = document.createElement('div');
        this.content.classList.add('widget-html-content');
        this.el.appendChild(this.content);
        this.update(); // Set defaults.
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        this.content.innerHTML = this.model.get('value');
        return super.update();
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do === 'focus') {
            this.content.focus();
        }
        else if (content.do === 'blur') {
            this.content.blur();
        }
    }
}
class HTMLMathModel extends StringModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'HTMLMathView', _model_name: 'HTMLMathModel' });
    }
}
class HTMLMathView extends StringView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-htmlmath');
        this.content = document.createElement('div');
        this.content.classList.add('widget-htmlmath-content');
        this.el.appendChild(this.content);
        this.update(); // Set defaults.
    }
    /**
     * Update the contents of this view
     */
    update() {
        this.content.innerHTML = this.model.get('value');
        this.typeset(this.content);
        return super.update();
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do === 'focus') {
            this.content.focus();
        }
        else if (content.do === 'blur') {
            this.content.blur();
        }
    }
}
class LabelModel extends StringModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'LabelView', _model_name: 'LabelModel' });
    }
}
class LabelView extends StringView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-label');
        this.update(); // Set defaults.
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        this.typeset(this.el, this.model.get('value'));
        return super.update();
    }
}
class TextareaModel extends StringModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'TextareaView', _model_name: 'TextareaModel', rows: null, continuous_update: true });
    }
}
class TextareaView extends StringView {
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-textarea');
        this.textbox = document.createElement('textarea');
        this.textbox.setAttribute('rows', '5');
        this.textbox.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.uuid)();
        this.textbox.classList.add('widget-input');
        this.el.appendChild(this.textbox);
        this.update(); // Set defaults.
        this.listenTo(this.model, 'change:placeholder', (model, value, options) => {
            this.update_placeholder(value);
        });
        this.update_placeholder();
        this.updateTooltip();
    }
    update_placeholder(value) {
        const v = value || this.model.get('placeholder');
        this.textbox.setAttribute('placeholder', v.toString());
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed.  The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update(options) {
        if (options === undefined || options.updated_view !== this) {
            this.textbox.value = this.model.get('value');
            let rows = this.model.get('rows');
            if (rows === null) {
                rows = '';
            }
            this.textbox.setAttribute('rows', rows);
            this.textbox.disabled = this.model.get('disabled');
        }
        this.updateTabindex();
        this.updateTooltip();
        return super.update();
    }
    updateTabindex() {
        if (!this.textbox) {
            return; // we might be constructing the parent
        }
        const tabbable = this.model.get('tabbable');
        if (tabbable === true) {
            this.textbox.setAttribute('tabIndex', '0');
        }
        else if (tabbable === false) {
            this.textbox.setAttribute('tabIndex', '-1');
        }
        else if (tabbable === null) {
            this.textbox.removeAttribute('tabIndex');
        }
    }
    updateTooltip() {
        if (!this.textbox)
            return; // we might be constructing the parent
        const title = this.model.get('tooltip');
        if (!title) {
            this.textbox.removeAttribute('title');
        }
        else if (this.model.get('description').length === 0) {
            this.textbox.setAttribute('title', title);
        }
    }
    events() {
        return {
            'keydown textarea': 'handleKeyDown',
            'keypress textarea': 'handleKeypress',
            'input textarea': 'handleChanging',
            'change textarea': 'handleChanged',
        };
    }
    /**
     * Handle key down
     *
     * Stop propagation so the event isn't sent to the application.
     */
    handleKeyDown(e) {
        e.stopPropagation();
    }
    /**
     * Handles key press
     *
     * Stop propagation so the keypress isn't sent to the application.
     */
    handleKeypress(e) {
        e.stopPropagation();
    }
    /**
     * Triggered on input change
     */
    handleChanging(e) {
        if (this.model.get('continuous_update')) {
            this.handleChanged(e);
        }
    }
    /**
     * Sync the value with the kernel.
     *
     * @param e Event
     */
    handleChanged(e) {
        const target = e.target;
        this.model.set('value', target.value, { updated_view: this });
        this.touch();
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do === 'focus') {
            this.textbox.focus();
        }
        else if (content.do === 'blur') {
            this.textbox.blur();
        }
    }
}
class TextModel extends StringModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'TextView', _model_name: 'TextModel', continuous_update: true });
    }
}
class TextView extends StringView {
    constructor() {
        super(...arguments);
        this.inputType = 'text';
    }
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('widget-text');
        this.textbox = document.createElement('input');
        this.textbox.setAttribute('type', this.inputType);
        this.textbox.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.uuid)();
        this.textbox.classList.add('widget-input');
        this.el.appendChild(this.textbox);
        this.update(); // Set defaults.
        this.listenTo(this.model, 'change:placeholder', (model, value, options) => {
            this.update_placeholder(value);
        });
        this.update_placeholder();
        this.updateTabindex();
        this.updateTooltip();
    }
    update_placeholder(value) {
        this.textbox.setAttribute('placeholder', value || this.model.get('placeholder'));
    }
    updateTabindex() {
        if (!this.textbox) {
            return; // we might be constructing the parent
        }
        const tabbable = this.model.get('tabbable');
        if (tabbable === true) {
            this.textbox.setAttribute('tabIndex', '0');
        }
        else if (tabbable === false) {
            this.textbox.setAttribute('tabIndex', '-1');
        }
        else if (tabbable === null) {
            this.textbox.removeAttribute('tabIndex');
        }
    }
    updateTooltip() {
        if (!this.textbox)
            return; // we might be constructing the parent
        const title = this.model.get('tooltip');
        if (!title) {
            this.textbox.removeAttribute('title');
        }
        else if (this.model.get('description').length === 0) {
            this.textbox.setAttribute('title', title);
        }
    }
    update(options) {
        /**
         * Update the contents of this view
         *
         * Called when the model is changed.  The model may have been
         * changed by another view or by a state update from the back-end.
         */
        if (options === undefined || options.updated_view !== this) {
            if (this.textbox.value !== this.model.get('value')) {
                this.textbox.value = this.model.get('value');
            }
            this.textbox.disabled = this.model.get('disabled');
        }
        return super.update();
    }
    events() {
        return {
            'keydown input': 'handleKeyDown',
            'keypress input': 'handleKeypress',
            'input input': 'handleChanging',
            'change input': 'handleChanged',
        };
    }
    /**
     * Handle key down
     *
     * Stop propagation so the keypress isn't sent to the application.
     */
    handleKeyDown(e) {
        e.stopPropagation();
    }
    /**
     * Handles text submission
     */
    handleKeypress(e) {
        e.stopPropagation();
        // The submit message is deprecated in widgets 7
        if (e.keyCode === 13) {
            // Return key
            this.send({ event: 'submit' });
        }
    }
    /**
     * Handles user input.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    handleChanging(e) {
        if (this.model.get('continuous_update')) {
            this.handleChanged(e);
        }
    }
    /**
     * Handles user input.
     *
     * Calling model.set will trigger all of the other views of the
     * model to update.
     */
    handleChanged(e) {
        const target = e.target;
        this.model.set('value', target.value, { updated_view: this });
        this.touch();
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do === 'focus') {
            this.textbox.focus();
        }
        else if (content.do === 'blur') {
            this.textbox.blur();
        }
    }
}
class PasswordModel extends TextModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'PasswordView', _model_name: 'PasswordModel' });
    }
}
class PasswordView extends TextView {
    constructor() {
        super(...arguments);
        this.inputType = 'password';
    }
}
/**
 * Combobox widget model class.
 */
class ComboboxModel extends TextModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'ComboboxModel', _view_name: 'ComboboxView', options: [], ensure_options: false });
    }
}
/**
 * Combobox widget view class.
 */
class ComboboxView extends TextView {
    constructor() {
        super(...arguments);
        this.isInitialRender = true;
    }
    render() {
        this.datalist = document.createElement('datalist');
        this.datalist.id = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.uuid)();
        super.render();
        this.textbox.setAttribute('list', this.datalist.id);
        this.el.appendChild(this.datalist);
        this.updateTooltip();
    }
    update(options) {
        super.update(options);
        if (!this.datalist) {
            return;
        }
        const valid = this.isValid(this.model.get('value'));
        this.highlightValidState(valid);
        // Check if we need to update options
        if ((options !== undefined && options.updated_view) ||
            (!this.model.hasChanged('options') && !this.isInitialRender)) {
            // Value update only, keep current options
            return;
        }
        this.isInitialRender = false;
        const opts = this.model.get('options');
        const optionFragment = document.createDocumentFragment();
        for (const v of opts) {
            const o = document.createElement('option');
            o.value = v;
            optionFragment.appendChild(o);
        }
        this.datalist.replaceChildren(...optionFragment.children);
    }
    isValid(value) {
        if (true === this.model.get('ensure_option')) {
            const options = this.model.get('options');
            if (options.indexOf(value) === -1) {
                return false;
            }
        }
        return true;
    }
    handleChanging(e) {
        // Override to validate value
        const target = e.target;
        const valid = this.isValid(target.value);
        this.highlightValidState(valid);
        if (valid) {
            super.handleChanging(e);
        }
    }
    handleChanged(e) {
        // Override to validate value
        const target = e.target;
        const valid = this.isValid(target.value);
        this.highlightValidState(valid);
        if (valid) {
            super.handleChanged(e);
        }
    }
    /**
     * Handle message sent to the front end.
     */
    handle_message(content) {
        if (content.do === 'focus') {
            this.textbox.focus();
        }
        else if (content.do === 'blur') {
            this.textbox.blur();
        }
    }
    highlightValidState(valid) {
        this.textbox.classList.toggle(INVALID_VALUE_CLASS, !valid);
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_tagsinput.js":
/*!*******************************************************!*\
  !*** ../../packages/controls/lib/widget_tagsinput.js ***!
  \*******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ColorsInputModel: () => (/* binding */ ColorsInputModel),
/* harmony export */   ColorsInputView: () => (/* binding */ ColorsInputView),
/* harmony export */   FloatsInputModel: () => (/* binding */ FloatsInputModel),
/* harmony export */   FloatsInputView: () => (/* binding */ FloatsInputView),
/* harmony export */   IntsInputModel: () => (/* binding */ IntsInputModel),
/* harmony export */   IntsInputView: () => (/* binding */ IntsInputView),
/* harmony export */   TagsInputModel: () => (/* binding */ TagsInputModel),
/* harmony export */   TagsInputView: () => (/* binding */ TagsInputView)
/* harmony export */ });
/* harmony import */ var d3_color__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! d3-color */ "../../node_modules/d3-color/src/color.js");
/* harmony import */ var d3_format__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! d3-format */ "../../node_modules/d3-format/src/defaultLocale.js");
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.




/**
 * Returns a new string after removing any leading and trailing whitespaces.
 * The original string is left unchanged.
 */
function trim(value) {
    return value.replace(/^\s+|\s+$/g, '');
}
/**
 * Clamp a number between min and max and return the result.
 */
function clamp(value, min, max) {
    return Math.min(Math.max(value, min), max);
}
/**
 * Remove children from an HTMLElement
 */
function removeChildren(el) {
    while (el.firstChild) {
        el.removeChild(el.firstChild);
    }
}
/**
 * Selection class which keeps track on selected indices.
 */
class Selection {
    constructor(start, dx, max) {
        this.start = start;
        this.dx = dx;
        this.max = max;
    }
    /**
     * Check if a given index is currently selected.
     */
    isSelected(index) {
        let min;
        let max;
        if (this.dx >= 0) {
            min = this.start;
            max = this.start + this.dx;
        }
        else {
            min = this.start + this.dx;
            max = this.start;
        }
        return min <= index && index < max;
    }
    /**
     * Update selection
     */
    updateSelection(dx) {
        this.dx += dx;
        if (this.start + this.dx > this.max) {
            this.dx = this.max - this.start;
        }
        if (this.start + this.dx < 0) {
            this.dx = -this.start;
        }
    }
}
class TagsInputBaseModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { value: [], placeholder: '\u200b', allowed_tags: null, allow_duplicates: true });
    }
}
class TagsInputBaseView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.DOMWidgetView {
    constructor() {
        super(...arguments);
        this.hoveredTag = null;
        this.hoveredTagIndex = null;
    }
    /**
     * Called when view is rendered.
     */
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('jupyter-widget-tagsinput');
        this.taginputWrapper = document.createElement('div');
        // The taginput is hidden until the user focuses on the widget
        // Unless there is no value
        if (this.model.get('value').length) {
            this.taginputWrapper.style.display = 'none';
        }
        else {
            this.taginputWrapper.style.display = 'inline-block';
        }
        this.datalistID = (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.uuid)();
        this.taginput = document.createElement('input');
        this.taginput.classList.add('jupyter-widget-tag');
        this.taginput.classList.add('jupyter-widget-taginput');
        this.taginput.setAttribute('list', this.datalistID);
        this.taginput.setAttribute('type', 'text');
        this.autocompleteList = document.createElement('datalist');
        this.autocompleteList.id = this.datalistID;
        this.updateAutocomplete();
        this.model.on('change:allowed_tags', this.updateAutocomplete.bind(this));
        this.updatePlaceholder();
        this.model.on('change:placeholder', this.updatePlaceholder.bind(this));
        this.taginputWrapper.classList.add('widget-text');
        this.taginputWrapper.appendChild(this.taginput);
        this.taginputWrapper.appendChild(this.autocompleteList);
        this.el.onclick = this.focus.bind(this);
        this.el.ondrop = (event) => {
            // Put the tag at the end of the list if there is no currently hovered tag
            const index = this.hoveredTagIndex == null ? this.tags.length : this.hoveredTagIndex;
            return this.ondrop(event, index);
        };
        this.el.ondragover = this.ondragover.bind(this);
        this.taginput.onchange = this.handleValueAdded.bind(this);
        this.taginput.oninput = this.resizeInput.bind(this);
        this.taginput.onkeydown = this.handleKeyEvent.bind(this);
        this.taginput.onblur = this.loseFocus.bind(this);
        this.resizeInput();
        this.inputIndex = this.model.get('value').length;
        this.selection = null;
        this.preventLoosingFocus = false;
        this.update();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update() {
        // Prevent hiding the input element and clearing the selection when updating everything
        this.preventLoosingFocus = true;
        removeChildren(this.el);
        this.tags = [];
        const value = this.model.get('value');
        this.inputIndex = value.length;
        for (const idx in value) {
            const index = parseInt(idx);
            const tag = this.createTag(value[index], index, this.selection != null && this.selection.isSelected(index));
            // Drag and drop
            tag.draggable = true;
            tag.ondragstart = ((index, value) => {
                return (event) => {
                    this.ondragstart(event, index, value, this.model.model_id);
                };
            })(index, value[index]);
            tag.ondrop = ((index) => {
                return (event) => {
                    this.ondrop(event, index);
                };
            })(index);
            tag.ondragover = this.ondragover.bind(this);
            tag.ondragenter = ((index) => {
                return (event) => {
                    this.ondragenter(event, index);
                };
            })(index);
            tag.ondragend = this.ondragend.bind(this);
            this.tags.push(tag);
            this.el.appendChild(tag);
        }
        this.el.insertBefore(this.taginputWrapper, this.el.children[this.inputIndex]);
        // The taginput is hidden until the user focuses on the widget
        // Unless there is no value
        if (this.model.get('value').length) {
            this.taginputWrapper.style.display = 'none';
        }
        else {
            this.taginputWrapper.style.display = 'inline-block';
        }
        this.preventLoosingFocus = false;
        return super.update();
    }
    /**
     * Update the auto-completion list
     */
    updateAutocomplete() {
        removeChildren(this.autocompleteList);
        const allowedTags = this.model.get('allowed_tags');
        for (const tag of allowedTags) {
            const option = document.createElement('option');
            option.value = tag;
            this.autocompleteList.appendChild(option);
        }
    }
    /**
     * Update the auto-completion list
     */
    updatePlaceholder() {
        this.taginput.placeholder = this.model.get('placeholder');
        this.resizeInput();
    }
    /**
     * Update the tags, called when the selection has changed and we need to update the tags CSS
     */
    updateTags() {
        const value = this.model.get('value');
        for (const idx in this.tags) {
            const index = parseInt(idx);
            this.updateTag(this.tags[index], value[index], index, this.selection != null && this.selection.isSelected(index));
        }
    }
    /**
     * Handle a new value is added from the input element
     */
    handleValueAdded(event) {
        const newTagValue = trim(this.taginput.value);
        const tagIndex = this.inputIndex;
        if (newTagValue == '') {
            return;
        }
        this.inputIndex++;
        const tagAdded = this.addTag(tagIndex, newTagValue);
        if (tagAdded) {
            // Clear the input and keep focus on it allowing the user to add more tags
            this.taginput.value = '';
            this.resizeInput();
            this.focus();
        }
    }
    /**
     * Add a new tag with a value of `tagValue` at the `index` position
     * Return true if the tag was correctly added, false otherwise
     */
    addTag(index, tagValue) {
        const value = this.model.get('value');
        let newTagValue;
        try {
            newTagValue = this.validateValue(tagValue);
        }
        catch (error) {
            return false;
        }
        const allowedTagValues = this.model.get('allowed_tags');
        if (allowedTagValues.length && !allowedTagValues.includes(newTagValue)) {
            // Do nothing for now, maybe show a proper error message?
            return false;
        }
        if (!this.model.get('allow_duplicates') && value.includes(newTagValue)) {
            // Do nothing for now, maybe add an animation to highlight the tag?
            return false;
        }
        // Clearing the current selection before setting the new value
        this.selection = null;
        // Making a copy so that backbone sees the change, and insert the new tag
        const newValue = [...value];
        newValue.splice(index, 0, newTagValue);
        this.model.set('value', newValue);
        this.model.save_changes();
        return true;
    }
    /**
     * Resize the input element
     */
    resizeInput() {
        let content;
        if (this.taginput.value.length != 0) {
            content = this.taginput.value;
        }
        else {
            content = this.model.get('placeholder');
        }
        const size = content.length + 1;
        this.taginput.setAttribute('size', String(size));
    }
    /**
     * Handle key events on the input element
     */
    handleKeyEvent(event) {
        const valueLength = this.model.get('value').length;
        // Do nothing if the user is typing something
        if (this.taginput.value.length) {
            return;
        }
        const currentElement = this.inputIndex;
        switch (event.key) {
            case 'ArrowLeft':
                if (event.ctrlKey && event.shiftKey) {
                    this.select(currentElement, -currentElement);
                }
                if (!event.ctrlKey && event.shiftKey) {
                    this.select(currentElement, -1);
                }
                if (event.ctrlKey) {
                    this.inputIndex = 0;
                }
                else {
                    this.inputIndex--;
                }
                break;
            case 'ArrowRight':
                if (event.ctrlKey && event.shiftKey) {
                    this.select(currentElement, valueLength - currentElement);
                }
                if (!event.ctrlKey && event.shiftKey) {
                    this.select(currentElement, 1);
                }
                if (event.ctrlKey) {
                    this.inputIndex = valueLength;
                }
                else {
                    this.inputIndex++;
                }
                break;
            case 'Backspace':
                if (this.selection) {
                    this.removeSelectedTags();
                }
                else {
                    this.removeTag(this.inputIndex - 1);
                }
                break;
            case 'Delete':
                if (this.selection) {
                    this.removeSelectedTags();
                }
                else {
                    this.removeTag(this.inputIndex);
                }
                break;
            default:
                // Do nothing by default
                return;
                break;
        }
        // Reset selection is shift key is not pressed
        if (!event.shiftKey) {
            this.selection = null;
        }
        this.inputIndex = clamp(this.inputIndex, 0, valueLength);
        this.update();
        this.focus();
    }
    /**
     * Function that gets called when a tag with a given `value` is being dragged.
     */
    ondragstart(event, index, tagValue, origin) {
        if (event.dataTransfer == null) {
            return;
        }
        event.dataTransfer.setData('index', String(index));
        event.dataTransfer.setData('tagValue', String(tagValue));
        event.dataTransfer.setData('origin', origin);
    }
    /**
     * Function that gets called when a tag has been dragged on the tag at the `index` position.
     */
    ondrop(event, index) {
        if (event.dataTransfer == null) {
            return;
        }
        event.preventDefault();
        event.stopPropagation();
        const draggedTagValue = event.dataTransfer.getData('tagValue');
        const draggedTagindex = parseInt(event.dataTransfer.getData('index'));
        const sameOrigin = event.dataTransfer.getData('origin') == this.model.model_id;
        // If something else than a tag was dropped, draggedTagindex should be NaN
        if (isNaN(draggedTagindex)) {
            return;
        }
        // If it's the same origin, the drag and drop results in a reordering
        if (sameOrigin) {
            const value = this.model.get('value');
            const newValue = [...value];
            // If the old position is on the left of the new position, we need to re-index the new position
            // after removing the tag at the old position
            if (draggedTagindex < index) {
                index--;
            }
            newValue.splice(draggedTagindex, 1); // Removing at the old position
            newValue.splice(index, 0, draggedTagValue); // Adding at the new one
            this.model.set('value', newValue);
            this.model.save_changes();
            return;
        }
        // Else we add a new tag with the given draggedTagValue
        this.addTag(index, draggedTagValue);
    }
    ondragover(event) {
        // This is needed for the drag and drop to work
        event.preventDefault();
    }
    ondragenter(event, index) {
        if (this.hoveredTag != null && this.hoveredTag != this.tags[index]) {
            this.hoveredTag.style.marginLeft = '1px';
        }
        this.hoveredTag = this.tags[index];
        this.hoveredTagIndex = index;
        this.hoveredTag.style.marginLeft = '30px';
    }
    ondragend() {
        if (this.hoveredTag != null) {
            this.hoveredTag.style.marginLeft = '1px';
        }
        this.hoveredTag = null;
        this.hoveredTagIndex = null;
    }
    /**
     * Select tags from `start` to `start + dx` not included.
     */
    select(start, dx) {
        const valueLength = this.model.get('value').length;
        if (!this.selection) {
            this.selection = new Selection(start, dx, valueLength);
        }
        else {
            this.selection.updateSelection(dx);
        }
    }
    /**
     * Remove all the selected tags.
     */
    removeSelectedTags() {
        const value = [...this.model.get('value')];
        const valueLength = value.length;
        // It is simpler to remove from right to left
        for (let idx = valueLength - 1; idx >= 0; idx--) {
            if (this.selection != null && this.selection.isSelected(idx)) {
                value.splice(idx, 1);
                // Move the input to the left if we remove a tag that is before the input
                if (idx < this.inputIndex) {
                    this.inputIndex--;
                }
            }
        }
        this.model.set('value', value);
        this.model.save_changes();
    }
    /**
     * Remove a tag given its index in the list
     */
    removeTag(tagIndex) {
        const value = [...this.model.get('value')];
        value.splice(tagIndex, 1);
        // Move the input to the left if we remove a tag that is before the input
        if (tagIndex < this.inputIndex) {
            this.inputIndex--;
        }
        this.model.set('value', value);
        this.model.save_changes();
    }
    /**
     * Focus on the input element
     */
    focus() {
        this.taginputWrapper.style.display = 'inline-block';
        this.taginput.focus();
    }
    /**
     * Lose focus on the input element
     */
    loseFocus() {
        if (this.preventLoosingFocus) {
            return;
        }
        // Only hide the input if we have tags displayed
        if (this.model.get('value').length) {
            this.taginputWrapper.style.display = 'none';
        }
        this.selection = null;
        this.updateTags();
    }
    preinitialize() {
        // Must set this before the initialize method creates the element
        this.tagName = 'div';
    }
    /**
     * Validate an input tag typed by the user, returning the correct tag type. This should be overridden in subclasses.
     */
    validateValue(value) {
        return value;
    }
}
class TagsInputModel extends TagsInputBaseModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { value: [], tag_style: '', _view_name: 'TagsInputView', _model_name: 'TagsInputModel' });
    }
}
class TagsInputView extends TagsInputBaseView {
    /**
     * Create the string tag
     */
    createTag(value, index, selected) {
        const tag = document.createElement('div');
        const style = this.model.get('tag_style');
        tag.classList.add('jupyter-widget-tag');
        tag.classList.add(TagsInputView.class_map[style]);
        if (selected) {
            tag.classList.add('mod-active');
        }
        tag.appendChild(document.createTextNode(this.getTagText(value)));
        const i = document.createElement('i');
        i.classList.add('fa');
        i.classList.add('fa-times');
        i.classList.add('jupyter-widget-tag-close');
        tag.appendChild(i);
        i.onmousedown = ((index) => {
            return () => {
                this.removeTag(index);
                this.loseFocus();
            };
        })(index);
        return tag;
    }
    /**
     * Returns the text that should be displayed in the tag element
     */
    getTagText(value) {
        return value;
    }
    /**
     * Update a given tag
     */
    updateTag(tag, value, index, selected) {
        if (selected) {
            tag.classList.add('mod-active');
        }
        else {
            tag.classList.remove('mod-active');
        }
    }
}
TagsInputView.class_map = {
    primary: 'mod-primary',
    success: 'mod-success',
    info: 'mod-info',
    warning: 'mod-warning',
    danger: 'mod-danger',
};
class ColorsInputModel extends TagsInputBaseModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { value: [], _view_name: 'ColorsInputView', _model_name: 'ColorsInputModel' });
    }
}
class ColorsInputView extends TagsInputBaseView {
    /**
     * Create the Color tag
     */
    createTag(value, index, selected) {
        const tag = document.createElement('div');
        const color = value;
        const darkerColor = d3_color__WEBPACK_IMPORTED_MODULE_2__["default"](value).darker().toString();
        tag.classList.add('jupyter-widget-tag');
        tag.classList.add('jupyter-widget-colortag');
        if (!selected) {
            tag.style.backgroundColor = color;
        }
        else {
            tag.classList.add('mod-active');
            tag.style.backgroundColor = darkerColor;
        }
        const i = document.createElement('i');
        i.classList.add('fa');
        i.classList.add('fa-times');
        i.classList.add('jupyter-widget-tag-close');
        tag.appendChild(i);
        i.onmousedown = ((index) => {
            return () => {
                this.removeTag(index);
                this.loseFocus();
            };
        })(index);
        return tag;
    }
    /**
     * Update a given tag
     */
    updateTag(tag, value, index, selected) {
        const color = value;
        const darkerColor = d3_color__WEBPACK_IMPORTED_MODULE_2__["default"](value).darker().toString();
        if (!selected) {
            tag.classList.remove('mod-active');
            tag.style.backgroundColor = color;
        }
        else {
            tag.classList.add('mod-active');
            tag.style.backgroundColor = darkerColor;
        }
    }
    /**
     * Validate an input tag typed by the user, returning the correct tag type. This should be overridden in subclasses.
     */
    validateValue(value) {
        if (d3_color__WEBPACK_IMPORTED_MODULE_2__["default"](value) == null) {
            throw value + ' is not a valid Color';
        }
        return value;
    }
}
class NumbersInputModel extends TagsInputModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { min: null, max: null });
    }
}
class NumbersInputView extends TagsInputView {
    render() {
        // Initialize text formatter
        this.model.on('change:format', () => {
            this.formatter = d3_format__WEBPACK_IMPORTED_MODULE_3__.format(this.model.get('format'));
            this.update();
        });
        this.formatter = d3_format__WEBPACK_IMPORTED_MODULE_3__.format(this.model.get('format'));
        super.render();
    }
    /**
     * Returns the text that should be displayed in the tag element
     */
    getTagText(value) {
        return this.formatter(this.parseNumber(value));
    }
    /**
     * Validate an input tag typed by the user, returning the correct tag type. This should be overridden in subclasses.
     */
    validateValue(value) {
        const parsed = this.parseNumber(value);
        const min = this.model.get('min');
        const max = this.model.get('max');
        if (isNaN(parsed) ||
            (min != null && parsed < min) ||
            (max != null && parsed > max)) {
            throw (value +
                ' is not a valid number, it should be in the range [' +
                min +
                ', ' +
                max +
                ']');
        }
        return parsed;
    }
}
class FloatsInputModel extends NumbersInputModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'FloatsInputView', _model_name: 'FloatsInputModel', format: '.1f' });
    }
}
class FloatsInputView extends NumbersInputView {
    parseNumber(value) {
        return parseFloat(value);
    }
}
class IntsInputModel extends NumbersInputModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _view_name: 'IntsInputView', _model_name: 'IntsInputModel', format: 'd' });
    }
}
class IntsInputView extends NumbersInputView {
    parseNumber(value) {
        const int = parseInt(value);
        if (int != parseFloat(value)) {
            throw value + ' should be an integer';
        }
        return int;
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_time.js":
/*!**************************************************!*\
  !*** ../../packages/controls/lib/widget_time.js ***!
  \**************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   TimeModel: () => (/* binding */ TimeModel),
/* harmony export */   TimeView: () => (/* binding */ TimeView),
/* harmony export */   deserialize_time: () => (/* binding */ deserialize_time),
/* harmony export */   serialize_time: () => (/* binding */ serialize_time),
/* harmony export */   time_serializers: () => (/* binding */ time_serializers)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "../../packages/controls/lib/utils.js");
/* harmony import */ var _widget_description__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_description */ "../../packages/controls/lib/widget_description.js");
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
// Copyright (c) Vidar Tonaas Fauske
// Distributed under the terms of the Modified BSD License.



const PARSER = /(\d\d):(\d\d)(:(\d\d)(.(\d{1,3})\d*)?)?/;
function serialize_time(value) {
    if (value === null) {
        return null;
    }
    else {
        const res = PARSER.exec(value);
        if (res === null) {
            return null;
        }
        return {
            hours: Math.min(23, parseInt(res[1], 10)),
            minutes: Math.min(59, parseInt(res[2], 10)),
            seconds: res[4] ? Math.min(59, parseInt(res[4], 10)) : 0,
            milliseconds: res[6] ? parseInt(res[6], 10) : 0,
        };
    }
}
function deserialize_time(value) {
    if (value === null) {
        return null;
    }
    else {
        const parts = [
            `${value.hours.toString().padStart(2, '0')}:${value.minutes
                .toString()
                .padStart(2, '0')}`,
        ];
        if (value.seconds > 0 || value.milliseconds > 0) {
            parts.push(`:${value.seconds.toString().padStart(2, '0')}`);
            if (value.milliseconds > 0) {
                parts.push(`.${value.milliseconds.toString().padStart(3, '0')}`);
            }
        }
        return parts.join('');
    }
}
const time_serializers = {
    serialize: serialize_time,
    deserialize: deserialize_time,
};
class TimeModel extends _widget_core__WEBPACK_IMPORTED_MODULE_2__.CoreDescriptionModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: TimeModel.model_name, _view_name: TimeModel.view_name, value: null, disabled: false, min: null, max: null, step: 60 });
    }
}
TimeModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_2__.CoreDescriptionModel.serializers), { value: time_serializers, min: time_serializers, max: time_serializers });
TimeModel.model_name = 'TimeModel';
TimeModel.view_name = 'TimeView';
class TimeView extends _widget_description__WEBPACK_IMPORTED_MODULE_1__.DescriptionView {
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-inline-hbox');
        this.el.classList.add('widget-timepicker');
        this._timepicker = document.createElement('input');
        this._timepicker.setAttribute('type', 'time');
        this._timepicker.id = this.label.htmlFor = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.uuid)();
        this.el.appendChild(this._timepicker);
        this.listenTo(this.model, 'change:value', this._update_value);
        this.listenTo(this.model, 'change', this.update2);
        this._update_value();
        this.update2();
    }
    /**
     * Update the contents of this view
     *
     * Called when the model is changed. The model may have been
     * changed by another view or by a state update from the back-end.
     */
    update2(model, options) {
        if (options === undefined || options.updated_view !== this) {
            this._timepicker.disabled = this.model.get('disabled');
            this._timepicker.min = this.model.get('min');
            this._timepicker.max = this.model.get('max');
            this._timepicker.step = this.model.get('step');
        }
        return super.update();
    }
    events() {
        // Typescript doesn't understand that these functions are called, so we
        // specifically use them here so it knows they are being used.
        void this._picker_change;
        void this._picker_focusout;
        return {
            'change [type="time"]': '_picker_change',
            'focusout [type="time"]': '_picker_focusout',
        };
    }
    _update_value(model, newValue, options) {
        if (options === undefined || options.updated_view !== this) {
            this._timepicker.value = this.model.get('value');
        }
    }
    _picker_change() {
        if (!this._timepicker.validity.badInput) {
            this.model.set('value', this._timepicker.value, { updated_view: this });
            this.touch();
        }
    }
    _picker_focusout() {
        if (this._timepicker.validity.badInput) {
            this.model.set('value', null, { updated_view: this });
            this.touch();
        }
    }
}


/***/ }),

/***/ "../../packages/controls/lib/widget_upload.js":
/*!****************************************************!*\
  !*** ../../packages/controls/lib/widget_upload.js ***!
  \****************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   FileUploadModel: () => (/* binding */ FileUploadModel),
/* harmony export */   FileUploadView: () => (/* binding */ FileUploadView)
/* harmony export */ });
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


class FileUploadModel extends _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'FileUploadModel', _view_name: 'FileUploadView', accept: '', description: 'Upload', disabled: false, icon: 'upload', button_style: '', multiple: false, value: [], error: '', style: null });
    }
}
FileUploadModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_0__.CoreDOMWidgetModel.serializers), { 
    // use a dummy serializer for value to circumvent the default serializer.
    value: { serialize: (x) => x } });
class FileUploadView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_1__.DOMWidgetView {
    preinitialize() {
        // Must set this before the initialize method creates the element
        this.tagName = 'button';
    }
    render() {
        super.render();
        this.el.classList.add('jupyter-widgets');
        this.el.classList.add('widget-upload');
        this.el.classList.add('jupyter-button');
        this.fileInput = document.createElement('input');
        this.fileInput.type = 'file';
        this.fileInput.style.display = 'none';
        this.el.addEventListener('click', () => {
            this.fileInput.click();
        });
        this.fileInput.addEventListener('click', () => {
            this.fileInput.value = '';
        });
        this.fileInput.addEventListener('change', () => {
            var _a;
            const promisesFile = [];
            Array.from((_a = this.fileInput.files) !== null && _a !== void 0 ? _a : []).forEach((file) => {
                promisesFile.push(new Promise((resolve, reject) => {
                    const fileReader = new FileReader();
                    fileReader.onload = () => {
                        // We know we can read the result as an array buffer since
                        // we use the `.readAsArrayBuffer` method
                        const content = fileReader.result;
                        resolve({
                            content,
                            name: file.name,
                            type: file.type,
                            size: file.size,
                            last_modified: file.lastModified,
                        });
                    };
                    fileReader.onerror = () => {
                        reject();
                    };
                    fileReader.onabort = fileReader.onerror;
                    fileReader.readAsArrayBuffer(file);
                }));
            });
            Promise.all(promisesFile)
                .then((files) => {
                this.model.set({
                    value: files,
                    error: '',
                });
                this.touch();
            })
                .catch((err) => {
                console.error('error in file upload: %o', err);
                this.model.set({
                    error: err,
                });
                this.touch();
            });
        });
        this.listenTo(this.model, 'change:button_style', this.update_button_style);
        this.set_button_style();
        this.update(); // Set defaults.
    }
    update() {
        this.el.disabled = this.model.get('disabled');
        this.el.setAttribute('title', this.model.get('tooltip'));
        const value = this.model.get('value');
        const description = `${this.model.get('description')} (${value.length})`;
        const icon = this.model.get('icon');
        if (description.length || icon.length) {
            this.el.textContent = '';
            if (icon.length) {
                const i = document.createElement('i');
                i.classList.add('fa');
                i.classList.add('fa-' + icon);
                if (description.length === 0) {
                    i.classList.add('center');
                }
                this.el.appendChild(i);
            }
            this.el.appendChild(document.createTextNode(description));
        }
        this.fileInput.accept = this.model.get('accept');
        this.fileInput.multiple = this.model.get('multiple');
        return super.update();
    }
    update_button_style() {
        this.update_mapped_classes(FileUploadView.class_map, 'button_style', this.el);
    }
    set_button_style() {
        this.set_mapped_classes(FileUploadView.class_map, 'button_style', this.el);
    }
}
FileUploadView.class_map = {
    primary: ['mod-primary'],
    success: ['mod-success'],
    info: ['mod-info'],
    warning: ['mod-warning'],
    danger: ['mod-danger'],
};


/***/ }),

/***/ "../../packages/controls/lib/widget_video.js":
/*!***************************************************!*\
  !*** ../../packages/controls/lib/widget_video.js ***!
  \***************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   VideoModel: () => (/* binding */ VideoModel),
/* harmony export */   VideoView: () => (/* binding */ VideoView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _widget_core__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget_core */ "../../packages/controls/lib/widget_core.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


class VideoModel extends _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'VideoModel', _view_name: 'VideoView', format: 'mp4', width: '', height: '', autoplay: true, loop: true, controls: true, value: new DataView(new ArrayBuffer(0)) });
    }
}
VideoModel.serializers = Object.assign(Object.assign({}, _widget_core__WEBPACK_IMPORTED_MODULE_1__.CoreDOMWidgetModel.serializers), { value: {
        serialize: (value) => {
            return new DataView(value.buffer.slice(0));
        },
    } });
class VideoView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    render() {
        /**
         * Called when view is rendered.
         */
        super.render();
        this.luminoWidget.addClass('jupyter-widgets');
        this.luminoWidget.addClass('widget-image');
        this.update(); // Set defaults.
    }
    update() {
        /**
         * Update the contents of this view
         *
         * Called when the model is changed.  The model may have been
         * changed by another view or by a state update from the back-end.
         */
        let url;
        const format = this.model.get('format');
        const value = this.model.get('value');
        if (format !== 'url') {
            const blob = new Blob([value], {
                type: `video/${this.model.get('format')}`,
            });
            url = URL.createObjectURL(blob);
        }
        else {
            url = new TextDecoder('utf-8').decode(value.buffer);
        }
        // Clean up the old objectURL
        const oldurl = this.el.src;
        this.el.src = url;
        if (oldurl) {
            URL.revokeObjectURL(oldurl);
        }
        // Height and width
        const width = this.model.get('width');
        if (width !== undefined && width.length > 0) {
            this.el.setAttribute('width', width);
        }
        else {
            this.el.removeAttribute('width');
        }
        const height = this.model.get('height');
        if (height !== undefined && height.length > 0) {
            this.el.setAttribute('height', height);
        }
        else {
            this.el.removeAttribute('height');
        }
        // Video attributes
        this.el.loop = this.model.get('loop');
        this.el.autoplay = this.model.get('autoplay');
        this.el.controls = this.model.get('controls');
        return super.update();
    }
    remove() {
        if (this.el.src) {
            URL.revokeObjectURL(this.el.src);
        }
        super.remove();
    }
    preinitialize() {
        // Must set this before the initialize method creates the element
        this.tagName = 'video';
    }
}


/***/ }),

/***/ "../../packages/controls/package.json":
/*!********************************************!*\
  !*** ../../packages/controls/package.json ***!
  \********************************************/
/***/ ((module) => {

module.exports = /*#__PURE__*/JSON.parse('{"name":"@jupyter-widgets/controls","version":"5.0.12","description":"Jupyter interactive widgets","repository":{"type":"git","url":"https://github.com/jupyter-widgets/ipywidgets.git"},"license":"BSD-3-Clause","author":"Project Jupyter","main":"lib/index.js","typings":"lib/index.d.ts","files":["lib/**/*.map","lib/**/*.d.ts","lib/**/*.js","css/*.css","dist/","src"],"scripts":{"build":"npm run build:src && npm run build:css","build:css":"lessc css/nouislider.less css/nouislider.css && postcss --use postcss-import --use postcss-cssnext -o css/widgets.built.css css/widgets.css","build:src":"tsc --build","build:test":"tsc --build test && webpack --config test/webpack.conf.js","clean":"npm run clean:src","clean:src":"rimraf lib && rimraf tsconfig.tsbuildinfo","prepublish":"npm run clean && npm run build","test":"npm run test:unit","test:coverage":"npm run build:test && webpack --config test/webpack-cov.conf.js && karma start test/karma-cov.conf.js","test:unit":"npm run test:unit:firefox && npm run test:unit:chrome","test:unit:chrome":"npm run test:unit:default -- --browsers=Chrome","test:unit:default":"npm run build:test && karma start test/karma.conf.js --log-level debug","test:unit:firefox":"npm run test:unit:default -- --browsers=Firefox","test:unit:firefox:headless":"npm run test:unit:default -- --browsers=FirefoxHeadless","test:unit:ie":"npm run test:unit:default -- --browsers=IE"},"dependencies":{"@jupyter-widgets/base":"^6.0.11","@lumino/algorithm":"^1 || ^2","@lumino/domutils":"^1 || ^2","@lumino/messaging":"^1 || ^2","@lumino/signaling":"^1 || ^2","@lumino/widgets":"^1 || ^2","d3-color":"^3.0.1","d3-format":"^3.0.1","jquery":"^3.1.1","nouislider":"15.4.0"},"devDependencies":{"@jupyterlab/services":"^6.0.0 || ^7.0.0","@types/d3-color":"^3.0.2","@types/d3-format":"^3.0.1","@types/expect.js":"^0.3.29","@types/jquery":"^3.5.16","@types/mathjax":"^0.0.37","@types/mocha":"^9.0.0","@types/node":"^17.0.2","chai":"^4.0.0","css-loader":"^6.5.1","expect.js":"^0.3.1","istanbul-instrumenter-loader":"^3.0.1","karma":"^6.3.3","karma-chrome-launcher":"^3.1.0","karma-coverage":"^2.0.3","karma-firefox-launcher":"^2.1.1","karma-ie-launcher":"^1.0.0","karma-mocha":"^2.0.1","karma-mocha-reporter":"^2.2.5","karma-webpack":"^5.0.0","less":"^4.1.2","mocha":"^9.0.0","npm-run-all":"^4.1.5","postcss":"^8.3.2","postcss-cli":"^9.1.0","postcss-cssnext":"^3.1.0","postcss-import":"^14.0.2","postcss-loader":"^6.1.0","rimraf":"^3.0.2","sinon":"^12.0.1","sinon-chai":"^3.3.0","style-loader":"^3.3.1","typescript":"~4.9.4","webpack":"^5.65.0"}}');

/***/ })

}]);
//# sourceMappingURL=packages_controls_lib_index_js.6d5d05e0ec5e0c9f8185.js.map